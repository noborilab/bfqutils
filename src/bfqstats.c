/*
 *   bfqstats: FastQ statistics
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

#define BFQSTATS_VERSION      "1.0"
#define BFQSTATS_YEAR        "2025"

#define DEFAULT_KMER              6
#define MAX_KMER                 10

#define error(do_exit, msg, ...) do { \
    fprintf(stderr, "[E::%s] " msg "\n", __func__, ##__VA_ARGS__); \
    if (do_exit) { \
      fputs("Encountered fatal error, exiting. Run bfqstats -h for usage.\n", stderr); \
      exit(EXIT_FAILURE); \
    } } while (0) 
#define quit(msg, ...) error(true, msg, ##__VA_ARGS__)

#define min(x, y) ((x) < (y) ? (x) : (y))

KSEQ_INIT(gzFile, gzread)

static void help(void) {
    printf(
        "bfqstats v%s  Copyright (C) %s  Benjamin Jean-Marie Tremblay\n"
        "\n"
        "Usage:  bfqstats [options] -i reads.fq[.gz]\n"
        " -i <file>  Input reads. Use '-' for stdin.\n"
        " -l <file>  Read Length histogram.\n"
        " -g <file>  Read GC content histogram.\n"
        " -q <file>  Mean read quality histogram.\n"
        " -Q <file>  Per-position mean quality histogram.\n"
        " -b <file>  Per-position base content histogram.\n"
        " -k <file>  K-mer counts and obs/exp ratios.\n"
        " -o <file>  Send summary stats to a file instead of stderr.\n"
        " -K <int>   K-mer size for -k. Default: %d\n"
        " -n <int>   Only examine this number of reads.\n"
        " -O         Send the reads to stdout.\n"
        " -z         If -O, compress as gzip.\n"
        " -N         Do not print summary stats to stderr.\n"
        " -v         Print the version and exit.\n"
        " -h         Print this help message and exit.\n"
        , BFQSTATS_VERSION, BFQSTATS_YEAR
        , DEFAULT_KMER
    );
}

static inline void *calloc_or_die(size_t size, const char *func_name) {
    void *result = calloc(1, size);
    if (result == NULL) {
        fprintf(stderr, "[E::%s] Out of memory (requested %lu B).\n", func_name, size);
        exit(EXIT_FAILURE);
    }
    return result;
}
#define alloc(size) calloc_or_die((size), __func__)

static const uint64_t pow4[16] = {
          1 ,          4 ,        16 ,         64
,       256 ,       1024 ,      4096 ,      16384
,     65536 ,     262144 ,   1048576 ,    4194304
,  16777216 ,   67108864 , 268435456 , 1073741824
};

static const unsigned char char2index[256] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 
};

typedef struct stat_t {
    int64_t  l;
    size_t   m;
    int64_t *d;
} stat_t;

typedef struct statf_t {
    int64_t  l;
    size_t   m;
    double  *d;
} statf_t;

static stat_t *initStat(const size_t m) {
    stat_t *stat = alloc(sizeof(stat_t));
    stat->d = alloc(sizeof(int64_t) * m);
    stat->m = m;
    return stat;
}

static statf_t *initStatf(const size_t m) {
    statf_t *stat = alloc(sizeof(statf_t));
    stat->d = alloc(sizeof(double) * m);
    stat->m = m;
    return stat;
}

static inline void allocMoreStat(stat_t *stat, const int l) {
  if ((l + 1) < stat->m) {
      stat->l = l;
  } else {
      size_t prev_m = stat->m;
      while (stat->m < (l + 1)) stat->m *= 2;
      stat->d = (int64_t *) realloc(stat->d, sizeof(int64_t) * stat->m);
      if (stat->d == NULL) quit("Out of memory (requested %lu B)", stat->m * sizeof(int64_t));
      stat->l = l;
      for (size_t i = prev_m; i < stat->m; i++) {
          stat->d[i] = 0;
      }
  }
}

static inline void allocMoreStatf(statf_t *stat, const int l) {
  if ((l + 1) < stat->m) {
      stat->l = l;
  } else {
      const size_t prev_m = stat->m;
      while (stat->m < (l + 1)) stat->m *= 2;
      stat->d = (double *) realloc(stat->d, sizeof(double) * stat->m);
      if (stat->d == NULL) quit("Out of memory (requested %lu B)", stat->m * sizeof(double));
      stat->l = l;
      for (size_t i = prev_m; i < stat->m; i++) {
          stat->d[i] = 0;
      }
  }
}

static FILE *newWriteFile(const char *f) {
    FILE *o = fopen(f, "w");
    if (o == NULL) quit("Error: Failed to create output file %s [%s]", f, strerror(errno));
    return o;
}

static inline int kmer2int(const unsigned char *seq, const int k, const int offset) {
    int kmer = 0;
    for (unsigned int j = 0, i = k - 1; i < -1; j++, i--) {
        const unsigned char ind = char2index[seq[offset + j]];
        if (ind == 4) return -1;
        kmer += pow4[i] * ind;
    }
    return kmer;
}

static void int2kmer(const int ind, const int k, char *kmer) {
    const char index2dna[5] = "ACGT";
    const int alphLen = 4;
    int decomposed[MAX_KMER];
    decomposed[k - 1] = ind;
    for (unsigned int i = k - 2; i < -1; i--) {
      decomposed[i] = decomposed[i + 1] / alphLen;
    }
    kmer[0] = index2dna[decomposed[0]];
    kmer[k] = '\0';
    for (int i = 1; i < k; i++) {
      kmer[i] = index2dna[alphLen + (decomposed[i] - (decomposed[i - 1] + 1) * alphLen)];
    }
}

static int statMedian(const stat_t *stat, const int64_t total) {
    const int64_t medVal = total / 2;
    int64_t rtotal = 0;
    for (int i = 0; i < stat->l; i++) {
        rtotal += stat->d[i] * (int64_t) (i + 1); 
        if (rtotal > medVal) return i + 1;
    }
    return -1;
}

typedef struct val_t {
    int v; int i;
} val_t;

static int cmpv(const void *a, const void *b) {
    const val_t *l = (val_t *) a;
    const val_t *r = (val_t *) b;
    if (l->v < r->v) return 1;
    if (l->v > r->v) return -1;
    return 0;
}

static int *statOrder(const stat_t *stat) {
    int *indices = alloc(stat->l * sizeof(int));
    val_t *v = alloc(stat->l * sizeof(val_t));
    for (int i = 0; i < stat->l; i++) {
        v[i].v = stat->d[i];
        v[i].i = i;
    }
    qsort(v, stat->l, sizeof(val_t), cmpv);
    for (int i = 0; i < stat->l; i++) {
        indices[i] = v[i].i;
    }
    free(v);
    return indices;
}

static int calcGC(const char *seq, const int l) {
    int gc = 0;
    for (int i = 0; i < l; i++) {
        gc += (seq[i] == 'G' | seq[i] == 'C');
    }
    return (int) round(100.0 * gc / l);
}

int main(int argc, char *argv[]) {
    
    setlocale(LC_NUMERIC, "en_US.UTF-8");  // For thousandths sep

    gzFile input = NULL;
    bool outReads = false, outReadsAreGz = false;
    gzFile outputGz = NULL;

    bool quiet = false, noStats = true;

    FILE *ovStats_f = NULL;
    FILE *lengthHist_f = NULL, *gcHist_f = NULL, *qualHist_f = NULL;
    FILE *ppQualHist_f = NULL, *ppBaseHist_f = NULL, *kCounts_f = NULL;
    int kmer = DEFAULT_KMER;
    int64_t nReadsMax = 0;

    char *end = NULL;

    int opt;
    while ((opt = getopt(argc, argv, "i:l:g:q:Q:b:k:o:K:n:OzNvh")) != -1) {
        switch (opt) {
          case 'i':
                if (optarg[0] == '-' && optarg[1] == '\0') {
                    input = gzdopen(fileno(stdin), "r");
                } else {
                    input = gzopen(optarg, "r");
                    if (input == NULL) quit("Failed to open file %s [%s]", optarg, strerror(errno));
                }
                break;
            case 'l':
                noStats = false;
                lengthHist_f = newWriteFile(optarg);
                break;
            case 'g':
                noStats = false;
                gcHist_f = newWriteFile(optarg);
                break;
            case 'q':
                noStats = false;
                qualHist_f = newWriteFile(optarg);
                break;
            case 'Q': 
                noStats = false;
                ppQualHist_f = newWriteFile(optarg);
                break;
            case 'b':
                noStats = false;
                ppBaseHist_f = newWriteFile(optarg);
                break;
            case 'k':
                noStats = false;
                kCounts_f = newWriteFile(optarg);
                break;
            case 'o':
                noStats = false;
                ovStats_f = newWriteFile(optarg);
                break;
            case 'K':
                kmer = atoi(optarg);
                if (kmer < 2 || kmer > MAX_KMER) {
                    quit("Error: -K must be an integer from 2 to %d", MAX_KMER);
                }
                break;
            case 'n':
                nReadsMax = (int64_t) strtod(optarg, &end);
                if (nReadsMax < 0) quit("Error: -n must be a positive integer");
                break;
            case 'O':
                outReads = true;
                break;
            case 'z':
                outReadsAreGz = true;
                outputGz = gzdopen(1, "wb");
                break;
            case 'N':
                quiet = true;
                break;
            case 'v':
                printf("bfqstats v%s\n", BFQSTATS_VERSION);
                exit(EXIT_SUCCESS);
            case 'h':
                help();
                exit(EXIT_SUCCESS);
            default:
                fputs("Encountered fatal error, exiting. Run bfqstats -h for usage.\n", stderr);
                exit(EXIT_FAILURE);
        }
    }

    if (ovStats_f == NULL) {
        ovStats_f = fopen("/dev/stdout", "w");
    }

    if (!quiet) noStats = false;

    if (input == NULL) {
        quit("Error: Missing -i");
    }

    if (!outReads && outReadsAreGz) {
        quit("Error: -z cannot be used in absence of -o");
    }

    stat_t *lengthHist = initStat(256);
    stat_t *gcHist = initStat(128);
    statf_t *ppQualHist = initStatf(256);

    stat_t *ppCountHist = initStat(256);

    stat_t *bAHist = initStat(256);
    stat_t *bCHist = initStat(256);
    stat_t *bGHist = initStat(256);
    stat_t *bTHist = initStat(256);
    stat_t *bNHist = initStat(256);

    stat_t *qualHist = initStat(128);
    stat_t *kCounts = initStat(pow4[kmer]);
    kCounts->l = pow4[kmer];

    kseq_t *read = kseq_init(input);
    int retVal;

    int64_t nReads = 0, bases = 0, q20 = 0, q30 = 0, nKmers = 0;
    int64_t nA = 0, nC = 0, nG = 0, nT = 0, nN = 0;
    double avgQual, baseQual;
    int readQual;

    while ((retVal = kseq_read(read)) >= 0) {
        if (nReadsMax && (nReads + 1) > nReadsMax) break;
        nReads++; 
        if (!noStats) {
            // TODO: try out read duplication
            if (read->seq.l + 1 > lengthHist->l) {
                allocMoreStat(lengthHist, read->seq.l + 1);
            }
            if (read->seq.l > ppQualHist->l) {
                allocMoreStatf(ppQualHist, read->seq.l);
                allocMoreStat(ppCountHist, read->seq.l);
                allocMoreStat(bAHist, read->seq.l);
                allocMoreStat(bCHist, read->seq.l);
                allocMoreStat(bGHist, read->seq.l);
                allocMoreStat(bTHist, read->seq.l);
                allocMoreStat(bNHist, read->seq.l);
            }
            lengthHist->d[read->seq.l]++;
            bases += read->seq.l;

            gcHist->d[calcGC(read->seq.s, read->seq.l)]++;

            avgQual = 0;
            for (int i = 0; i < read->seq.l; i++) {
                if ((read->qual.s[i] - 33) > 19) {
                    q20++;
                    if ((read->qual.s[i] - 33) > 29) {
                        q30++;
                    }
                }
                baseQual = pow(10, ((double) (read->qual.s[i] - 33) / (-10.0)));
                avgQual += baseQual;
                ppQualHist->d[i] += baseQual;
                ppCountHist->d[i]++;
                switch (read->seq.s[i]) {
                    case 'a':
                    case 'A':
                        bAHist->d[i]++; nA++; break;
                    case 'c':
                    case 'C':
                        bCHist->d[i]++; nC++; break;
                    case 'g':
                    case 'G':
                        bGHist->d[i]++; nG++; break;
                    case 't':
                    case 'T':
                        bTHist->d[i]++; nT++; break;
                    default:
                        bNHist->d[i]++; nN++; break;
                }
                if (1 + (read->seq.l - i) > kmer) {
                    const int ind = kmer2int((unsigned char *) read->seq.s, kmer, i);
                    if (ind != -1) {
                        nKmers++;
                        kCounts->d[ind]++;
                    }
                }
            }
            avgQual /= (double) read->seq.l;
            readQual = min(94, (int) (log10(avgQual) * -10.0));
            if (readQual + 1 > qualHist->l) allocMoreStat(qualHist, readQual + 1);
            qualHist->d[readQual]++;
        }
        if (outReads) {
            if (outReadsAreGz) {
                gzprintf(outputGz, "@%s %s\n%s\n+\n%s\n", read->name.s, read->comment.s,
                    read->seq.s, read->qual.s);
            } else {
                printf("@%s %s\n%s\n+\n%s\n", read->name.s, read->comment.s,
                    read->seq.s, read->qual.s);
            }
        }
    }

    kseq_destroy(read);

    gzclose(input);
    if (outReadsAreGz) gzclose(outputGz);
    if (!nReads) quit("Error: Found 0 reads in input");

    if (!noStats) {

        const int lengthMedian = statMedian(lengthHist, bases);
        statf_t *kEnr = initStatf(kCounts->m);
        kEnr->l = kCounts->l;
        const double kExp = (double) nKmers / kCounts->l;
        for (int i = 0; i < kCounts->l; i++) {
            kEnr->d[i] = (double) kCounts->d[i] / kExp;
        }
        int *kOrd = statOrder(kCounts);

        fprintf(ovStats_f, "Total reads: %'lld\n", nReads);
        fprintf(ovStats_f, "Median length: %d\n", lengthMedian);
        fprintf(ovStats_f, "Base composition:\n A: %.1f%%\n C: %.1f%%\n G: %.1f%%\n T: %.1f%%\n N: %.1f%%\n",
            100.0 * (double) nA / bases, 100.0 * (double) nC / bases,
            100.0 * (double) nG / bases, 100.0 * (double) nT / bases,
            100.0 * (double) nN / bases);
        fprintf(ovStats_f, "Q20 bases: %.1f%%\n", 100.0 * (double) q20 / bases);
        fprintf(ovStats_f, "Q30 bases: %.1f%%\n", 100.0 * (double) q30 / bases);
        fprintf(ovStats_f, "Top 10 k-mers:\n");

        char kmerStr[MAX_KMER + 1];
        for (int i = 0; i < 10; i++) {
            int2kmer(kOrd[i], kmer, kmerStr);
            fprintf(ovStats_f, " %s: %'lld (obs/exp = %.1f)\n", kmerStr, kCounts->d[kOrd[i]], kEnr->d[kOrd[i]]);
        }

        if (kCounts_f != NULL) {
            fprintf(kCounts_f, "Kmer\tCount\tObsExp\n");
            for (int i = 0; i < kCounts->l; i++) {
                int2kmer(i, kmer, kmerStr);
                fprintf(kCounts_f, "%s\t%lld\t%f\n", kmerStr, kCounts->d[i], kEnr->d[i]);
            }
        }
        free(kEnr->d); free(kEnr); free(kOrd);

        if (lengthHist_f != NULL) {
            fprintf(lengthHist_f, "Length\tCount\tFraction\n");
            for (int i = 0; i < lengthHist->l; i++) {
                fprintf(lengthHist_f, "%d\t%lld\t%f\n", i, lengthHist->d[i], (double) lengthHist->d[i] / nReads);
            }
        }
        if (gcHist_f != NULL) {
            fprintf(gcHist_f, "PctGC\tCount\tFraction\n");
            for (int i = 0; i < 101; i++) {
                fprintf(gcHist_f, "%d%%\t%lld\t%f\n", i, gcHist->d[i], (double) gcHist->d[i] / nReads);
            }
        }
        if (qualHist_f != NULL) {
            fprintf(qualHist_f, "MeanQ\tCount\tFraction\n");
            for (int i = 0; i < qualHist->l; i++) {
                fprintf(qualHist_f, "%d\t%lld\t%f\n", i, qualHist->d[i], (double) qualHist->d[i] / nReads);
            }
        }
        if (ppQualHist_f != NULL) {
            fprintf(ppQualHist_f, "Pos\tMeanQ\n");
            double qual;
            for (int i = 0; i < ppQualHist->l; i++) {
                qual = ppQualHist->d[i];
                qual /= (double) ppCountHist->d[i];
                fprintf(ppQualHist_f, "%d\t%d\n", i + 1, (int) (log10(qual) * -10.0));
            }
        }
        if (ppBaseHist_f != NULL) {
            fprintf(ppBaseHist_f, "Pos\tA\tC\tG\tT\tN\n");
            for (int i = 0; i < bAHist->l; i++) {
                fprintf(ppBaseHist_f, "%d\t%lld\t%lld\t%lld\t%lld\t%lld\n", i + 1,
                    bAHist->d[i], bCHist->d[i], bGHist->d[i], bTHist->d[i], bNHist->d[i]);
            }
        }

    }

    fclose(ovStats_f);
    fclose(lengthHist_f);
    fclose(gcHist_f);
    fclose(qualHist_f);
    fclose(ppQualHist_f);
    fclose(ppBaseHist_f);
    fclose(kCounts_f);

    free(lengthHist->d); free(lengthHist);
    free(kCounts->d); free(kCounts);
    free(qualHist->d); free(qualHist);
    free(ppQualHist->d); free(ppQualHist); free(ppCountHist);
    free(bAHist->d); free(bAHist);
    free(bCHist->d); free(bCHist);
    free(bGHist->d); free(bGHist);
    free(bTHist->d); free(bTHist);
    free(bNHist->d); free(bNHist);

    return EXIT_SUCCESS;

}

