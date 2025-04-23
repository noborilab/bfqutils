# bfqutils: Ben's FastQ Utilities

* FastQ PE read merging: [bfqmerge](#bfqmerge)
* FastQ trimming for single-end reads: [bfqtrimse](#bfqtrimse)
* FastQ statistics: [bfqstats](#bfqstats)

## Installation

Requires gcc/clang and GNU Make, tested with macOS and Linux. bfqutils is written in C and comes bundled with the Cloudflare fork of Zlib.

```sh
git clone https://github.com/noborilab/bfqutils
cd bfqutils
make libz 
make release
make test     # optional
```

To dynamically link to a system Zlib instead of building a static library:

```sh
make z_dyn=1 release
```

## bfqmerge

This tool is meant to be used as a drop-in replacement for [`fastp --merge`](https://github.com/OpenGene/fastp?tab=readme-ov-file#merge-paired-end-reads). The only differences from fastp are a smaller default minimum merge length (since this tool is primarily intended for libraries with small inserts) and more aggressive base correction. Apart from that, bfqmerge is intended to be a lightweight replacement with extremely low memory usage (<2 MB) and very good single-thread performance (<15 min on 2x100M PE150 reads, ~10 min without the `-z` flag). Remember to change the default values in `src/bfqmerge.c` if that would better suit your primary use case.

### Usage

```
bfqmerge v1.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay

Usage:  bfqmerge [options] R1.fq[.gz] R2.fq[.gz] > merged.fq
 -o <int>  Required overlap for a merge to occur. Default: 15
 -d <int>  Maximum number of mismatches between alignments. Default: 5
 -p <dbl>  Maximum fraction of mismatches between alignments. Default: 0.2
 -Q <int>  Minimum PHRED+33 quality to consider a base high quality. Default: 15
 -u <dbl>  Maximum fraction of bases allowed to be low quality. Default: 0.4
 -n <int>  Maximum number of Ns allowed. Default: 5
 -g <int>  Number of Gs to trigger polyG tail trimming. Default: 10
 -t <int>  Mean window quality threshold for trimming 3-prime bases. Default: 20.
 -w <int>  Window size of 3-prime base trimming. Default: 5
 -m <int>  Max merged read length.
 -z        Compress the output as gzip.
 -q        Make the program quiet.
 -v        Print the version and exit.
 -h        Print this help message and exit.
 ```

### Implementation details

bfqmerge is a simple program, which performs the following operations:

1. Detect and trim low quality bases from the 3' end of the reads (`-t`, `-w`). This helps reduce the possibility of false positive alignments between incorrect read segments, especially when expecting short inserts.

2. Detect and trim polyG tails (`-g`). Failure to remove these can lead to bfqmerge thinking that sufficiently long stretches of Gs found in both reads of a pair are the overlapping parts of the reads. Due to the way some machines work, sequencing past the end of a read can lead to long stretches of Gs with high quality scores (which won't be trimmed in step 1).

3. Find the best overlap, depending on user settings (`-o`, `-d`, `-p`).

4. If no overlap is found, discard the reads. Otherwise, create a new merged read. For each overlapping position, the base and quality score is taken from whichever of the two reads has a better score in that position.

5. Check if the merged read passes quality filters (`-Q`, `-u`, `-n`). If yes, then compress the output if desired (`-z`), then write to `stdout`.

With the exception of `-o`, all default values are identical to fastp. These generally work quite well, though be aware that it is impossible to get read merging right 100% of the time. False positive merging events and false negative read discards will almost always occur for any combination of settings given a sufficiently diverse set of reads.

## bfqtrimse

A simple FastQ trimming and quality filtering tool for single-end reads. The order of operations is similar to bfqmerge, replacing overlap analysis with adapter sequence matching. Use '-' for stdin.

### Usage

```
bfqtrimse v1.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay

Usage:  bfqtrimse [options] reads.fq[.gz] > trimmed.fq
 -a <str>   Adapter sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
 -Q <int>   Minimum PHRED+33 quality to consider a base high quality. Default: 15
 -u <dbl>   Maximum fraction of bases allowed to be low quality. Default: 0.4
 -n <int>   Maximum number of Ns allowed. Default: 5
 -g <int>   Number of Gs to trigger polyG tail trimming. Default: 10
 -t <int>   Mean window quality threshold for trimming 3-prime bases. Default: 20.
 -w <int>   Window size of 3-prime base trimming. Default: 5
 -M <int>   Min trimmed read length. Default: 15
 -m <int>   Max trimmed read length.
 -z         Compress the output as gzip.
 -q         Make the program quiet.
 -v         Print the version and exit.
 -h         Print this help message and exit.
```

## bfqstats

Meant to be used in a pipe with bfqmerge or bfqtrimse. The default top enriched K-mer summary can help spot bad trimming/merging. Use '-' for stdin.

### Usage

```
bfqstats v1.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay

Usage:  bfqstats [options] reads.fq[.gz]
 -l <file>  Read Length histogram.
 -g <file>  Read GC content histogram.
 -q <file>  Mean read quality histogram.
 -Q <file>  Per-position mean quality histogram.
 -b <file>  Per-position base content histogram.
 -k <file>  K-mer counts and obs/exp ratios.
 -o <file>  Send summary stats to a file instead of stderr.
 -K <int>   K-mer size for -k. Default: 6
 -n <int>   Only examine this number of reads.
 -O         Send the reads to stdout.
 -z         If -O, compress as gzip.
 -N         Do not print summary stats to stderr.
 -v         Print the version and exit.
 -h         Print this help message and exit.
```

