/*
 *   bfqtrimse: FastQ adapter trimming for single-end reads
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

#define BFQTRIMSE_VERSION             "1.0"
#define BFQTRIMSE_YEAR               "2025"

#define DEFAULT_MIN_PHRED_QUAL          15
#define DEFAULT_MAX_NON_QUALIFIED      0.4
#define DEFAULT_MAX_N_BASES              5
#define DEFAULT_POLYG_N                 10
#define DEFAULT_TRIM_QUAL               20
#define DEFAULT_TRIM_WIN                 5
#define DEFAULT_MIN_LEN                 15

#define DEFAULT_ADP "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

#define error(do_exit, msg, ...) do { \
    fprintf(stderr, "[E::%s] " msg "\n", __func__, ##__VA_ARGS__); \
    if (do_exit) { \
      fputs("Encountered fatal error, exiting. Run bfqtrimse -h for usage.\n", stderr); \
      exit(EXIT_FAILURE); \
    } } while (0) 
#define quit(msg, ...) error(true, msg, ##__VA_ARGS__)

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

KSEQ_INIT(gzFile, gzread)

static void help(void) {
    printf(
        "bfqtrimse v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "\n"
        "Usage:  bfqtrimse [options] R1.fq[.gz] > trimmed.fq \n"
        " -a <str>   Adapter sequence for SE trimming. Default: %s\n"
        " -Q <int>   Minimum PHRED+33 quality to consider a base high quality. Default: %d\n"
        " -u <dbl>   Maximum fraction of bases allowed to be low quality. Default: %.1f\n"
        " -n <int>   Maximum number of Ns allowed. Default: %d\n"
        " -g <int>   Number of Gs to trigger polyG tail trimming. Default: %d\n"
        " -t <int>   Mean window quality threshold for trimming 3-prime bases. Default: %d.\n"
        " -w <int>   Window size of 3-prime base trimming. Default: %d\n"
        " -M <int>   Min trimmed read length. Default: %d\n"
        " -m <int>   Max trimmed read length.\n"
        " -z         Compress the output as gzip.\n"
        " -q         Make the program quiet.\n"
        " -v         Print the version and exit.\n"
        " -h         Print this help message and exit.\n"
        , BFQTRIMSE_VERSION, BFQTRIMSE_YEAR
        , DEFAULT_ADP
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

static inline int trim_read(kseq_t *read, const char *adp, const int alen, const int cstart) {

    // https://github.com/OpenGene/fastp/blob/master/src/adaptertrimmer.cpp

    const int allowOneMismatchForEach = 8;
    const int matchReq = 4;

    int pos = 0;
    bool found = false;
    int start = cstart;

    for (pos = start; pos < (int) read->seq.l - matchReq; pos++) {
        int cmplen = min((int) read->seq.l - pos, alen);
        int allowedMismatch = cmplen / allowOneMismatchForEach;
        int mismatch = 0;
        bool matched = true;
        for (int i = max(0, -pos); i < cmplen; i++) {
            if (adp[i] != read->seq.s[i + pos]) {
                mismatch++;
                if (mismatch > allowedMismatch) {
                    matched = false;
                    break;
                }
            }
        }
        if (matched) {
            found = true;
            break;
        }
    }

    if (found) {
        if (pos < 0) {
            return 0;
        } else {
            read->seq.s[pos] = 0;
            read->seq.l = (size_t) pos;
            read->qual.s[pos] = 0;
            read->qual.l = (size_t) pos;
        }
        return pos;
    }

    return 0;

}

int main(int argc, char *argv[]) {

    setlocale(LC_NUMERIC, "en_US.UTF-8");  // For thousandths sep

    int minPhredQual = DEFAULT_MIN_PHRED_QUAL + 33;
    float maxNonQualified = DEFAULT_MAX_NON_QUALIFIED; 
    int maxNBases = DEFAULT_MAX_N_BASES; 
    int polyG_n = DEFAULT_POLYG_N;
    int trimQ = DEFAULT_TRIM_QUAL;
    int trimW = DEFAULT_TRIM_WIN;
    int maxLen = 0, minLen = 0;
    char *adp = DEFAULT_ADP;
    bool gzip = false, quiet = false;
    
    int opt;
    while ((opt = getopt(argc, argv, "a:Q:u:n:t:w:m:M:zqvh")) != -1) {
        switch (opt) {
            case 'a':
                adp = optarg;
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
                if (minLen < 1) {
                    quit("Error: -M must be a positive integer");
                }
                break;
            case 'z':
                gzip = true;
                break;
            case 'q':
                quiet = true;
                break;
            case 'v':
                printf("bfqtrimse v%s\n", BFQTRIMSE_VERSION);
                exit(EXIT_SUCCESS);
            case 'h':
                help();
                exit(EXIT_SUCCESS);
            default:
                fputs("Encountered fatal error, exiting. Run bfqtrimse -h for usage.\n", stderr);
                exit(EXIT_FAILURE);
        }
    }

    const int n_files = argc - optind;
    if (!n_files) quit("Missing input files.");

    if (n_files != 1) {
        quit("Expected one input file, found %d.\n", n_files);
    }

    const int alen = strlen(adp);
    int cstart = 0;
    if (alen >= 16) {
        cstart = -4;
    } else if (alen >= 12) {
        cstart = -3;
    } else if (alen >= 8) {
        cstart = -2;
    }
    if (alen < 4) {
        quit("Error: Adapter sequence should be at least 4 bases long.");
    }

    int64_t nReads = 0, nTrimmed = 0;

    gzFile read_f;
    if (argv[optind][0] == '-' && argv[optind][1] == '\0') {
        read_f = gzdopen(fileno(stdin), "r");
    } else {
        read_f = gzopen(argv[optind], "r");
        if (read_f == NULL) quit("Failed to open file %s [%s]", argv[optind], strerror(errno));
    }

    gzFile gz;
    if (gzip) gz = gzdopen(1, "wb");
    /* gzFile read_f = gzopen(argv[optind], "r");  */
    kseq_t *read = kseq_init(read_f);
    int retVal;

    while ((retVal = kseq_read(read)) >= 0) {
        nReads++;
        if (!trim3p(read, trimQ, trimW)) continue;
        if (!trimPolyG(read, polyG_n)) continue;
        if (!trim_read(read, adp, alen, cstart)) continue;
        int lowQualNum = 0, nBaseNum = 0;
        for (int i = 0; i < (int) read->qual.l; i++) {
          if (read->qual.s[i] < minPhredQual) lowQualNum++;
          if (dnatable[(unsigned char) read->seq.s[i]] == 110) nBaseNum++;
        }
        if ((lowQualNum > (int) (maxNonQualified * (float) read->seq.l)) ||
            (nBaseNum > maxNBases) || (read->seq.l < minLen) ||
            (maxLen > 0 && read->seq.l > maxLen)) {
            continue;
        } else if (gzip) {
            gzprintf(gz, "@%s %s\n%s\n+\n%s\n", read->name.s,
                read->comment.l ? read->comment.s : "",
                read->seq.s, read->qual.s);
        } else {
            printf("@%s %s\n%s\n+\n%s\n", read->name.s,
                read->comment.l ? read->comment.s : "",
                read->seq.s, read->qual.s);
        }
        nTrimmed++;
    }

    kseq_destroy(read);
    gzclose(read_f);
    if (gzip) gzclose(gz);

    if (!quiet) {
        fprintf(stderr, "Processed %'lld reads\n", nReads);
        fprintf(stderr, "Output %'lld trimmed reads (%.2f%%)\n",
            nTrimmed, 100 * (double) nTrimmed / (double) nReads);
    }

    return EXIT_SUCCESS;

}

