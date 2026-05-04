/*
 *   bfqutils trimpe: FastQ PE adapter trimming via overlap analysis
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
#include "version.h"

#define DEFAULT_OVERLAP_REQUIRE         15
#define DEFAULT_DIFF_LIMIT               5
#define DEFAULT_DIFF_PERCENT_LIMIT     0.2
#define DEFAULT_MIN_PHRED_QUAL          15
#define DEFAULT_MAX_NON_QUALIFIED      0.4
#define DEFAULT_MAX_N_BASES              5
#define DEFAULT_POLYG_N                 10
#define DEFAULT_TRIM_QUAL               20
#define DEFAULT_TRIM_WIN                 5
#define DEFAULT_MIN_LEN                 15

#define error(do_exit, msg, ...) do { \
    fprintf(stderr, "[E::%s] " msg "\n", __func__, ##__VA_ARGS__); \
    if (do_exit) { \
      fputs("Encountered fatal error, exiting. Run bfqutils trimpe -h for usage.\n", stderr); \
      exit(EXIT_FAILURE); \
    } } while (0)
#define quit(msg, ...) error(true, msg, ##__VA_ARGS__)

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

KSEQ_INIT(gzFile, gzread)

static void help(void) {
    printf(
        "bfqutils v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "Usage:  bfqutils trimpe [options] R1.fq[.gz] R2.fq[.gz]\n"
        " -o <int>   Required overlap to detect adapters. Default: %d\n"
        " -d <int>   Maximum number of mismatches in the overlap. Default: %d\n"
        " -p <dbl>   Maximum fraction of mismatches in the overlap. Default: %.1f\n"
        " -1 <str>   Output file for trimmed R1 reads.\n"
        " -2 <str>   Output file for trimmed R2 reads.\n"
        " -i         Interleaved output to stdout. Mutually exclusive with -1/-2.\n"
        " -s <str>   Output file for singleton reads (surviving mate of a failed pair).\n"
        " -Q <int>   Minimum PHRED+33 quality to consider a base high quality. Default: %d\n"
        " -u <dbl>   Maximum fraction of bases allowed to be low quality. Default: %.1f\n"
        " -n <int>   Maximum number of Ns allowed. Default: %d\n"
        " -g <int>   Number of Gs to trigger polyG tail trimming. Default: %d\n"
        " -t <int>   Mean window quality threshold for trimming 3-prime bases. Default: %d.\n"
        " -w <int>   Window size of 3-prime base trimming. Default: %d\n"
        " -m <int>   Max read length.\n"
        " -M <int>   Min read length. Default: %d\n"
        " -z         Compress the output as gzip.\n"
        " -q         Make the program quiet.\n"
        " -v         Print the version and exit.\n"
        " -h         Print this help message and exit.\n"
        , BFQUTILS_VERSION, BFQUTILS_YEAR
        , DEFAULT_OVERLAP_REQUIRE
        , DEFAULT_DIFF_LIMIT
        , DEFAULT_DIFF_PERCENT_LIMIT
        , DEFAULT_MIN_PHRED_QUAL
        , DEFAULT_MAX_NON_QUALIFIED
        , DEFAULT_MAX_N_BASES
        , DEFAULT_POLYG_N
        , DEFAULT_TRIM_QUAL
        , DEFAULT_TRIM_WIN
        , DEFAULT_MIN_LEN
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

static int str_split(char *str, char **res, const char sep, const int max_size) {
  const char delim[2] = { sep, '\0' };
  int field_n = 0;
  if (strlen(str) == 0) return 0;
  res[field_n++] = strtok(str, delim);
  while (field_n < max_size && res[field_n - 1] != NULL) {
    res[field_n++] = strtok(NULL, delim);
  }
  if (field_n == max_size && res[field_n - 1] != NULL) return -1;
  return field_n - 1;
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
            if ((log10(avgQual) * -10.0) < (double) trimQ) {
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

static void fqtrimpe(const char *fwd, const char *rev,
    const int overlapRequire, const int diffLimit, const float diffPercentLimit,
    const bool gzip, const bool quiet, const char minPhredQual,
    const float maxNonQualified, const int maxNBases, const int polyG_n,
    const int trimQ, const int trimW, const int maxLen, const int minLen,
    gzFile gz_r1, gzFile gz_r2, gzFile gz_s,
    FILE *fp_r1,  FILE *fp_r2,  FILE *fp_s) {

    const int complete_compare_require = 50;
    int offset, overlap_len, diff;

    gzFile fwd_f = gzopen(fwd, "r");
    gzFile rev_f = gzopen(rev, "r");
    if (fwd_f == NULL) quit("Failed to open file %s [%s]", fwd, strerror(errno));
    if (rev_f == NULL) quit("Failed to open file %s [%s]", rev, strerror(errno));

    kseq_t *fwd_kseq = kseq_init(fwd_f);
    kseq_t *rev_kseq = kseq_init(rev_f);

    int fwd_ret_val, rev_ret_val;
    int64_t reads_n = 0, out_n = 0, sing_n = 0;

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
            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)) {
                goto finish_overlap;
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
            if (diff <= overlapDiffLimit || (diff > overlapDiffLimit && i > complete_compare_require)) {
                goto finish_overlap;
            }
            offset -= 1;
        }

        goto filter;

finish_overlap:;
        {
            const int from_fwd = overlap_len + max(0, offset);
            const int from_rev = offset > 0 ? (int) rev_kseq->seq.l - overlap_len : 0;
            const int insert_len = from_fwd + from_rev;
            if ((int) fwd_kseq->seq.l > insert_len) {
                fwd_kseq->seq.s[insert_len] = 0;  fwd_kseq->seq.l  = insert_len;
                fwd_kseq->qual.s[insert_len] = 0; fwd_kseq->qual.l = insert_len;
            }
            if ((int) rev_kseq->seq.l > insert_len) {
                rev_kseq->seq.s[insert_len] = 0;  rev_kseq->seq.l  = insert_len;
                rev_kseq->qual.s[insert_len] = 0; rev_kseq->qual.l = insert_len;
            }
        }

filter:;
        // https://github.com/OpenGene/fastp/blob/master/src/filter.cpp
        int fwd_lowQ = 0, fwd_nBase = 0;
        for (int i = 0; i < (int) fwd_kseq->qual.l; i++) {
            if (fwd_kseq->qual.s[i] < minPhredQual) fwd_lowQ++;
            if (fwd_kseq->seq.s[i] == 'N') fwd_nBase++;
        }
        const bool fwd_ok = (fwd_lowQ <= (int)(maxNonQualified * (float) fwd_kseq->seq.l))
            && (fwd_nBase <= maxNBases)
            && ((int) fwd_kseq->seq.l >= minLen)
            && (maxLen == 0 || (int) fwd_kseq->seq.l <= maxLen);

        int rev_lowQ = 0, rev_nBase = 0;
        for (int i = 0; i < (int) rev_kseq->qual.l; i++) {
            if (rev_kseq->qual.s[i] < minPhredQual) rev_lowQ++;
            if (rev_kseq->seq.s[i] == 'N') rev_nBase++;
        }
        const bool rev_ok = (rev_lowQ <= (int)(maxNonQualified * (float) rev_kseq->seq.l))
            && (rev_nBase <= maxNBases)
            && ((int) rev_kseq->seq.l >= minLen)
            && (maxLen == 0 || (int) rev_kseq->seq.l <= maxLen);

        if (fwd_ok && rev_ok) {
            out_n++;
            if (gzip) {
                gzprintf(gz_r1, "@%s %s\n%s\n+\n%s\n", fwd_kseq->name.s,
                    fwd_kseq->comment.l ? fwd_kseq->comment.s : "",
                    fwd_kseq->seq.s, fwd_kseq->qual.s);
                gzprintf(gz_r2, "@%s %s\n%s\n+\n%s\n", rev_kseq->name.s,
                    rev_kseq->comment.l ? rev_kseq->comment.s : "",
                    rev_kseq->seq.s, rev_kseq->qual.s);
            } else {
                fprintf(fp_r1, "@%s %s\n%s\n+\n%s\n", fwd_kseq->name.s,
                    fwd_kseq->comment.l ? fwd_kseq->comment.s : "",
                    fwd_kseq->seq.s, fwd_kseq->qual.s);
                fprintf(fp_r2, "@%s %s\n%s\n+\n%s\n", rev_kseq->name.s,
                    rev_kseq->comment.l ? rev_kseq->comment.s : "",
                    rev_kseq->seq.s, rev_kseq->qual.s);
            }
        } else if (fwd_ok && (gz_s || fp_s)) {
            sing_n++;
            if (gzip) {
                gzprintf(gz_s, "@%s %s\n%s\n+\n%s\n", fwd_kseq->name.s,
                    fwd_kseq->comment.l ? fwd_kseq->comment.s : "",
                    fwd_kseq->seq.s, fwd_kseq->qual.s);
            } else {
                fprintf(fp_s, "@%s %s\n%s\n+\n%s\n", fwd_kseq->name.s,
                    fwd_kseq->comment.l ? fwd_kseq->comment.s : "",
                    fwd_kseq->seq.s, fwd_kseq->qual.s);
            }
        } else if (rev_ok && (gz_s || fp_s)) {
            sing_n++;
            if (gzip) {
                gzprintf(gz_s, "@%s %s\n%s\n+\n%s\n", rev_kseq->name.s,
                    rev_kseq->comment.l ? rev_kseq->comment.s : "",
                    rev_kseq->seq.s, rev_kseq->qual.s);
            } else {
                fprintf(fp_s, "@%s %s\n%s\n+\n%s\n", rev_kseq->name.s,
                    rev_kseq->comment.l ? rev_kseq->comment.s : "",
                    rev_kseq->seq.s, rev_kseq->qual.s);
            }
        }
    }

    kseq_destroy(fwd_kseq);
    kseq_destroy(rev_kseq);

    gzclose(fwd_f);
    gzclose(rev_f);

    if (!quiet) {
        fprintf(stderr, "Processed %'lld read pairs\n", reads_n);
        fprintf(stderr, "Output %'lld read pairs (%.2f%%)\n",
            out_n, reads_n > 0 ? 100 * (double) out_n / (double) reads_n : 0.0);
        if (gz_s || fp_s) {
            fprintf(stderr, "Output %'lld singletons (%.2f%%)\n",
                sing_n, reads_n > 0 ? 100 * (double) sing_n / (double) reads_n : 0.0);
        }
    }

}

int main_trimpe(int argc, char *argv[]) {

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
    int maxLen = 0, minLen = DEFAULT_MIN_LEN;
    char *out1 = NULL, *out2 = NULL, *outs = NULL;
    bool interleaved = false;
    bool gzip = false, quiet = false;

    int opt;
    while ((opt = getopt(argc, argv, "o:d:p:1:2:s:iQ:u:n:g:t:w:m:M:zqvh")) != -1) {
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
            case '1':
                out1 = optarg;
                break;
            case '2':
                out2 = optarg;
                break;
            case 's':
                outs = optarg;
                break;
            case 'i':
                interleaved = true;
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
            case 'M':
                minLen = atoi(optarg);
                if (minLen < 0) {
                    quit("Error: -M must be greater or equal to 0");
                }
                break;
            case 'z':
                gzip = true;
                break;
            case 'q':
                quiet = true;
                break;
            case 'v':
                printf("bfqutils v%s\n", BFQUTILS_VERSION);
                exit(EXIT_SUCCESS);
            case 'h':
                help();
                exit(EXIT_SUCCESS);
            default:
                fputs("Encountered fatal error, exiting. Run bfqutils trimpe -h for usage.\n", stderr);
                exit(EXIT_FAILURE);
        }
    }

    const int n_files = argc - optind;
    if (!n_files) quit("Missing input files.");

    if (n_files != 2) {
        quit("Expected two input files, found %d.\n", n_files);
    }

    if (interleaved && (out1 || out2)) {
        quit("Error: -i is mutually exclusive with -1 / -2");
    }
    if (!interleaved && (!out1 || !out2)) {
        quit("Error: must specify -1 and -2, or -i");
    }

    char *fwd_f[1024];
    char *rev_f[1024];
    const int n_fwd = str_split(argv[optind], fwd_f, ',', 1024);
    const int n_rev = str_split(argv[optind + 1], rev_f, ',', 1024);

    if (n_fwd != n_rev) {
        fprintf(stderr, "Error: found %d R1 files and %d R2 files.\n", n_fwd, n_rev);
        exit(EXIT_FAILURE);
    }

    gzFile gz_r1 = NULL, gz_r2 = NULL, gz_s = NULL;
    FILE  *fp_r1 = NULL, *fp_r2 = NULL, *fp_s = NULL;

    if (interleaved) {
        if (gzip) {
            gz_r1 = gz_r2 = gzdopen(1, "wb");
        } else {
            fp_r1 = fp_r2 = stdout;
        }
    } else {
        if (gzip) {
            gz_r1 = gzopen(out1, "wb");
            gz_r2 = gzopen(out2, "wb");
            if (!gz_r1) quit("Failed to open output file %s [%s]", out1, strerror(errno));
            if (!gz_r2) quit("Failed to open output file %s [%s]", out2, strerror(errno));
        } else {
            fp_r1 = fopen(out1, "w");
            fp_r2 = fopen(out2, "w");
            if (!fp_r1) quit("Failed to open output file %s [%s]", out1, strerror(errno));
            if (!fp_r2) quit("Failed to open output file %s [%s]", out2, strerror(errno));
        }
    }

    if (outs) {
        if (gzip) {
            gz_s = gzopen(outs, "wb");
            if (!gz_s) quit("Failed to open singletons file %s [%s]", outs, strerror(errno));
        } else {
            fp_s = fopen(outs, "w");
            if (!fp_s) quit("Failed to open singletons file %s [%s]", outs, strerror(errno));
        }
    }

    for (int i = 0; i < n_fwd; i++) {
        fqtrimpe(fwd_f[i], rev_f[i], overlapRequire, diffLimit, diffPercentLimit,
            gzip, quiet, (char) minPhredQual, maxNonQualified, maxNBases, polyG_n,
            trimQ, trimW, maxLen, minLen,
            gz_r1, gz_r2, gz_s, fp_r1, fp_r2, fp_s);
    }

    if (gzip) {
        gzclose(gz_r1);
        if (!interleaved) gzclose(gz_r2);
        if (gz_s) gzclose(gz_s);
    } else {
        if (!interleaved) {
            fclose(fp_r1);
            fclose(fp_r2);
        }
        if (fp_s) fclose(fp_s);
    }

    return EXIT_SUCCESS;

}
