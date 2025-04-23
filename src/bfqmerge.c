/*
 *   bfqmerge: FastQ PE read merging
 *   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <locale.h>
#include <limits.h>
#include <errno.h>
#include <math.h>
#include "kseq.h"
#include "zlib.h"

#define BFQMERGE_VERSION             "1.0"
#define BFQMERGE_YEAR               "2025"

#define DEFAULT_OVERLAP_REQUIRE         15
#define DEFAULT_DIFF_LIMIT               5
#define DEFAULT_DIFF_PERCENT_LIMIT     0.2
#define DEFAULT_MIN_PHRED_QUAL          15
#define DEFAULT_MAX_NON_QUALIFIED      0.4
#define DEFAULT_MAX_N_BASES              5
#define DEFAULT_POLYG_N                 10
#define DEFAULT_TRIM_QUAL               20
#define DEFAULT_TRIM_WIN                 5

#define error(do_exit, msg, ...) do { \
    fprintf(stderr, "[E::%s] " msg "\n", __func__, ##__VA_ARGS__); \
    if (do_exit) { \
      fputs("Encountered fatal error, exiting. Run bfqmerge -h for usage.\n", stderr); \
      exit(EXIT_FAILURE); \
    } } while (0) 
#define quit(msg, ...) error(true, msg, ##__VA_ARGS__)

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

KSEQ_INIT(gzFile, gzread)

static inline void *calloc_or_die(size_t size, const char *func_name) {
    void *result = calloc(1, size);
    if (result == NULL) {
        fprintf(stderr, "[E::%s] Out of memory (requested %lu B).\n", func_name, size);
        exit(EXIT_FAILURE);
    }
    return result;
}
#define alloc(size) calloc_or_die((size), __func__)

static void help(void) {
    printf(
        "bfqmerge v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "\n"
        "Usage:  bfqmerge [options] R1.fq[.gz] R2.fq[.gz] > merged.fq\n"
        " -o <int>   Required overlap for a merge to occur. Default: %d\n"
        " -d <int>   Maximum number of mismatches between alignments. Default: %d\n"
        " -p <dbl>   Maximum fraction of mismatches between alignments. Default: %.1f\n"
        " -Q <int>   Minimum PHRED+33 quality to consider a base high quality. Default: %d\n"
        " -u <dbl>   Maximum fraction of bases allowed to be low quality. Default: %.1f\n"
        " -n <int>   Maximum number of Ns allowed. Default: %d\n"
        " -g <int>   Number of Gs to trigger polyG tail trimming. Default: %d\n"
        " -t <int>   Mean window quality threshold for trimming 3-prime bases. Default: %d.\n"
        " -w <int>   Window size of 3-prime base trimming. Default: %d\n"
        " -m <int>   Max merged read length.\n"
        " -z         Compress the output as gzip.\n"
        " -q         Make the program quiet.\n"
        " -v         Print the version and exit.\n"
        " -h         Print this help message and exit.\n"
        , BFQMERGE_VERSION, BFQMERGE_YEAR
        , DEFAULT_OVERLAP_REQUIRE
        , DEFAULT_DIFF_LIMIT
        , DEFAULT_DIFF_PERCENT_LIMIT
        , DEFAULT_MIN_PHRED_QUAL
        , DEFAULT_MAX_NON_QUALIFIED
        , DEFAULT_MAX_N_BASES
        , DEFAULT_POLYG_N
        , DEFAULT_TRIM_QUAL
        , DEFAULT_TRIM_WIN
    );
}

static const unsigned char dnatable[256] = {
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110,  65, 110,  67, 110, 110, 110,  71, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110,  84, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110,  65, 110,  67, 110, 110, 110,  71, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110,  84, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110,
    110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110
};

static const unsigned char rctable[256] = {
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78,
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78
};

static kstring_t *kstring_init(void) {
    kstring_t *read = alloc(sizeof(kstring_t));
    read->m = 256;
    read->s = alloc(sizeof(char) * read->m);
    return read;
}

static inline void merge_reads(const kseq_t *fwd_kseq, const kseq_t *rev_kseq, kseq_t *merged_kseq, const int offset, const int overlap_len) {

    const int from_fwd = overlap_len + max(0, offset);
    const int from_rev = offset > 0 ? (int) rev_kseq->seq.l - overlap_len : 0;
    merged_kseq->seq.l = from_fwd + from_rev;
    merged_kseq->qual.l = merged_kseq->seq.l;

    if ((merged_kseq->seq.l + 1) >= merged_kseq->seq.m) {
        while (merged_kseq->seq.m < (merged_kseq->seq.l + 1)) merged_kseq->seq.m *= 2;
        merged_kseq->qual.m = merged_kseq->seq.m;
        merged_kseq->seq.s = (char *) realloc(merged_kseq->seq.s, merged_kseq->seq.m);
        merged_kseq->qual.s = (char *) realloc(merged_kseq->qual.s, merged_kseq->qual.m);
    }

    for (int i = 0; i < from_fwd; i++) {
        merged_kseq->seq.s[i] = fwd_kseq->seq.s[i];    
        merged_kseq->qual.s[i] = fwd_kseq->qual.s[i];    
    }
    for (int i = 0; i < from_rev; i++) {
        merged_kseq->seq.s[from_fwd + i] = rctable[(unsigned char) rev_kseq->seq.s[(from_rev - i) - 1]]; 
        merged_kseq->qual.s[from_fwd + i] = rev_kseq->qual.s[(from_rev - i) - 1];
    }
    merged_kseq->seq.s[merged_kseq->seq.l] = 0;
    merged_kseq->qual.s[merged_kseq->qual.l] = 0;

    merged_kseq->name.l = fwd_kseq->name.l;
    if (merged_kseq->name.m < (merged_kseq->name.l + 1)) {
        while (merged_kseq->name.m < (merged_kseq->name.l + 1)) merged_kseq->name.m *= 2;
        merged_kseq->name.s = (char *) realloc(merged_kseq->name.s, merged_kseq->name.m);
    }
    for (int i = 0; i < (int) fwd_kseq->name.l; i++) {
        merged_kseq->name.s[i] = fwd_kseq->name.s[i];
    } 
    merged_kseq->name.s[fwd_kseq->name.l] = 0;

    merged_kseq->comment.l = fwd_kseq->comment.l;
    if (merged_kseq->comment.m < (merged_kseq->comment.l + 1)) {
        while (merged_kseq->comment.m < (merged_kseq->comment.l + 1)) merged_kseq->comment.m *= 2;
        merged_kseq->comment.s = (char *) realloc(merged_kseq->comment.s, merged_kseq->comment.m);
    }
    for (int i = 0; i < (int) fwd_kseq->comment.l; i++) {
        merged_kseq->comment.s[i] = fwd_kseq->comment.s[i];
    } 
    merged_kseq->comment.s[fwd_kseq->comment.l] = 0;

    const int start1 = max(0, offset);
    const int start2 = rev_kseq->seq.l - max(0, -offset) - 1;
    for (int i = 0; i < overlap_len; i++) {
        const int p1 = start1 + i;
        const int p2 = start2 - i;
        if ((unsigned char) rev_kseq->qual.s[p2] > (unsigned char) fwd_kseq->qual.s[p1]) {
            merged_kseq->qual.s[p1] = rev_kseq->qual.s[p2]; 
            merged_kseq->seq.s[p1] = rctable[(unsigned char) rev_kseq->seq.s[p2]];
        }
    }

}

static inline int trimPolyG(kseq_t *read, const int polyG_n) {
    // https://github.com/OpenGene/fastp/blob/master/src/polyx.cpp
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;
    int mismatch = 0, firstGPos = (int) read->seq.l - 1;
    int i, allowedMismatch;
    for (i = 0; i < (int) read->seq.l; i++) {
        if (read->seq.s[(int) read->seq.l - i - 1] != 'G') {
            mismatch++;
        } else {
            firstGPos = (int) read->seq.l - i - 1;
        }
        allowedMismatch = (i + 1) / allowOneMismatchForEach;
        if (mismatch > maxMismatch || (mismatch > allowedMismatch && i >= polyG_n - 1)) {
            break;
        }
    }
    if (i >= polyG_n) {
       read->seq.s[firstGPos] = 0; 
       read->seq.l = firstGPos + 1;
       read->qual.s[firstGPos] = 0; 
       read->qual.l = firstGPos + 1;
    }
    return (int) read->seq.l;
}

static inline int trim3p(kseq_t *read, const int trimQ, const int trimW) {
    double avgQual = 0;
    int toTrim = -1;
    if (read->seq.l >= trimW) {
        for (int i = read->seq.l - trimW; i > -1; i--) {
            avgQual = 0;
            for (int j = i; j < i + trimW; j++) {
                avgQual += pow(10, ((double) (read->qual.s[j] - 33) / (-10.0)));
            }
            avgQual /= trimW;
            if (((int) log10(avgQual) * -10.0) < trimQ) {
                toTrim = i;
            } else {
                break;
            }
        }
        if (toTrim > -1) {
            read->seq.s[toTrim] = 0;
            read->qual.s[toTrim] = 0;
            read->seq.l = toTrim;
            read->qual.l = toTrim;
        }
    }
    return (int) read->seq.l;
}

static void fqmerge(const char *fwd, const char *rev, const int overlapRequire, const int diffLimit, const float diffPercentLimit, const bool gzip, const bool quiet, const char minPhredQual, const float maxNonQualified, const int maxNBases, const int polyG_n, const int trimQ, const int trimW, const int maxLen) {

    gzFile gz;
    if (gzip) gz = gzdopen(1, "wb");

    const int complete_compare_require = 50;  // overlapDiffLimit only applies in first 50 bp?
    int offset, overlap_len, diff;

    gzFile fwd_f = gzopen(fwd, "r"); 
    gzFile rev_f = gzopen(rev, "r"); 
    if (fwd_f == NULL) quit("Failed to open file %s [%s]", fwd, strerror(errno));
    if (rev_f == NULL) quit("Failed to open file %s [%s]", rev, strerror(errno));

    kseq_t *fwd_kseq = kseq_init(fwd_f);
    kseq_t *rev_kseq = kseq_init(rev_f);
    kseq_t *merged_kseq = alloc(sizeof(kseq_t));
    merged_kseq->name = *kstring_init();
    merged_kseq->comment = *kstring_init();
    merged_kseq->seq = *kstring_init();
    merged_kseq->qual = *kstring_init();

    int fwd_ret_val, rev_ret_val;
    int64_t reads_n = 0, merged_n = 0, filtered_n = 0;

    while ((fwd_ret_val = kseq_read(fwd_kseq)) >= 0) {
        reads_n++;
        if ((rev_ret_val = kseq_read(rev_kseq)) < 0) {
            quit("PE FastQ files do not have the same number of reads.\n");
        }

        if (!trim3p(fwd_kseq, trimQ, trimW)) continue;
        if (!trim3p(rev_kseq, trimQ, trimW)) continue;
        if (!trimPolyG(fwd_kseq, polyG_n)) continue;
        if (!trimPolyG(rev_kseq, polyG_n)) continue;

        // https://github.com/OpenGene/fastp/blob/master/src/overlapanalysis.cpp
        offset = 0;
        while (offset <= (int) fwd_kseq->seq.l - overlapRequire) {
            overlap_len = min((int) fwd_kseq->seq.l - offset, (int) rev_kseq->seq.l);
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));
            diff = 0;
            int i;
            for (i = 0; i < overlap_len; i++) {
                const unsigned char fwd_c = dnatable[(unsigned char) fwd_kseq->seq.s[offset + i]];
                const unsigned char rev_c = rctable[(unsigned char) rev_kseq->seq.s[(int) rev_kseq->seq.l - i - 1]];
                if (fwd_c != rev_c) {
                    diff += 1;
                    if (diff > overlapDiffLimit && i < complete_compare_require) break;
                }
            }
            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)){
                merge_reads(fwd_kseq, rev_kseq, merged_kseq, offset, overlap_len);
                goto finish_merge;
            }
            offset += 1;
        }

        offset = 0;
        while (offset >= -((int) rev_kseq->seq.l - overlapRequire)) {
            overlap_len = min((int) fwd_kseq->seq.l, (int) rev_kseq->seq.l - abs(offset));
            int overlapDiffLimit = min(diffLimit, (int)(overlap_len * diffPercentLimit));
            diff = 0;
            int i;
            for (i = 0; i < overlap_len; i++) {
                const unsigned char fwd_c = dnatable[(unsigned char) fwd_kseq->seq.s[i]];
                const unsigned char rev_c = rctable[(unsigned char) rev_kseq->seq.s[(int) rev_kseq->seq.l - (-offset + i) - 1]];
                if (fwd_c != rev_c) {
                    diff += 1;
                    if (diff > overlapDiffLimit && i < complete_compare_require) break;
                }
            }
            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)){
                merge_reads(fwd_kseq, rev_kseq, merged_kseq, offset, overlap_len);
                goto finish_merge;
            }
            offset -= 1;
        }

        continue;

finish_merge:;
        const int from_fwd = overlap_len + max(0, offset);
        const int from_rev = offset > 0 ? (int) rev_kseq->seq.l - overlap_len : 0;
        // https://github.com/OpenGene/fastp/blob/master/src/filter.cpp
        int lowQualNum = 0, nBaseNum = 0;
        for (int i = 0; i < merged_kseq->qual.l; i++) {
          if (merged_kseq->qual.s[i] < minPhredQual) lowQualNum++;
          if (merged_kseq->seq.s[i] == 'N') nBaseNum++;
        }
        if ((lowQualNum > (int) (maxNonQualified * (float) merged_kseq->seq.l)) ||
            (nBaseNum > maxNBases) || (merged_kseq->seq.l < overlapRequire) ||
            (maxLen > 0 && merged_kseq->seq.l > maxLen)) {
            filtered_n++;
        } else if (gzip) {
            gzprintf(gz, "@%s %s merged_%d_%d\n%s\n+\n%s\n", merged_kseq->name.s, merged_kseq->comment.s,
                from_fwd, from_rev, merged_kseq->seq.s, merged_kseq->qual.s);
        } else {
            printf("@%s %s merged_%d_%d\n%s\n+\n%s\n", merged_kseq->name.s, merged_kseq->comment.s,
                from_fwd, from_rev, merged_kseq->seq.s, merged_kseq->qual.s);
        }
        merged_n++;
    }

    kseq_destroy(fwd_kseq);
    kseq_destroy(rev_kseq);
    kseq_destroy(merged_kseq);

    gzclose(fwd_f);
    gzclose(rev_f);

    if (gzip) gzclose(gz);

    if (!quiet) {
        fprintf(stderr, "Processed %'lld read pairs\n", reads_n);
        fprintf(stderr, "Created %'lld merged reads (%.2f%%)\n",
            merged_n, 100 * (double) merged_n / (double) reads_n);
        fprintf(stderr, "Output %'lld filtered reads (%.2f%%)\n",
            merged_n - filtered_n, 100 * (double) (merged_n - filtered_n) / (double) reads_n);
    }

}

int main(int argc, char *argv[]) {

    setlocale(LC_NUMERIC, "en_US.UTF-8");  // For thousandths sep

    int overlapRequire = DEFAULT_OVERLAP_REQUIRE;
    int diffLimit = DEFAULT_DIFF_LIMIT;
    float diffPercentLimit = DEFAULT_DIFF_PERCENT_LIMIT;
    int minPhredQual = DEFAULT_MIN_PHRED_QUAL + 33;
    float maxNonQualified = DEFAULT_MAX_NON_QUALIFIED; 
    int maxNBases = DEFAULT_MAX_N_BASES; 
    int polyG_n = DEFAULT_POLYG_N;
    int trimQ = DEFAULT_TRIM_QUAL;
    int trimW = DEFAULT_TRIM_WIN;
    int maxLen = 0;
    bool gzip = false, quiet = false;

    int opt;
    while ((opt = getopt(argc, argv, "o:d:p:Q:u:n:t:w:m:zqvh")) != -1) {
        switch (opt) {
            case 'o':
                overlapRequire = atoi(optarg);
                if (overlapRequire < 1) quit("Error: -o must be a positive integer");
                break;
            case 'd':
                diffLimit = atoi(optarg);
                if (diffLimit < 0) quit("Error: -d must be greater or equal to 0");
                break;
            case 'p':
                diffPercentLimit = atof(optarg);
                if (diffPercentLimit < 0.0f || diffPercentLimit > 1.0f) {
                    quit("Error: -p must be between 0 and 1");
                }
                break;
            case 'Q':
                minPhredQual = atoi(optarg);
                if (minPhredQual > 127 - 33) {
                    minPhredQual = 127 - 33;
                }
                if (minPhredQual < 0) {
                    minPhredQual = 0;
                }
                minPhredQual += 33;
                break;
            case 'u':
                maxNonQualified = atof(optarg);
                if (maxNonQualified < 0.0f || maxNonQualified > 1.0f) {
                    quit("Error: -u must be between 0 and 1");
                }
                break;
            case 'n':
                maxNBases = atoi(optarg);
                if (maxNBases < 0) {
                    quit("Error: -n must be greater or equal to 0");
                }
                break;
            case 'g':
                polyG_n = atoi(optarg);
                if (polyG_n < 0) {
                    quit("Error: -g must be greater or equal to 0");
                }
                break;
            case 't':
                trimQ = atoi(optarg);
                if (trimQ > 127 - 33) {
                    trimQ = 127 - 33;
                }
                if (trimQ < 0) {
                    trimQ = 0;
                }
                break;
            case 'w':
                trimW = atoi(optarg);
                if (trimW < 1) {
                    quit("Error: -w must be a positive integer");
                }
                break;
            case 'm':
                maxLen = atoi(optarg);
                if (maxLen < 1) {
                    quit("Error: -m must be a positive integer");
                }
                break;
            case 'z':
                gzip = true;
                break;
            case 'q':
                quiet = true;
                break;
            case 'v':
                printf("bfqmerge v%s\n", BFQMERGE_VERSION);
                exit(EXIT_SUCCESS);
            case 'h':
                help();
                exit(EXIT_SUCCESS);
            default:
                fputs("Encountered fatal error, exiting. Run bfqmerge -h for usage.\n", stderr);
                exit(EXIT_FAILURE);
        }
    }

    const int n_files = argc - optind;
    if (!n_files) quit("Missing input files.");

    if (n_files != 2) {
        quit("Expected two input files, found %d.\n", n_files);
    }

    fqmerge(argv[optind], argv[optind + 1], overlapRequire, diffLimit, diffPercentLimit, gzip, quiet, (char) minPhredQual, maxNonQualified, maxNBases, polyG_n, trimQ, trimW, maxLen);

    return EXIT_SUCCESS;

}

