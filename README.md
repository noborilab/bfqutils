# bfqutils: Ben's FastQ Utilities

* FastQ merging for PE reads: [bfqmerge](#bfqmerge)
* FastQ trimming for SE reads: [bfqtrimse](#bfqtrimse)

## Installation

Requires gcc/clang and GNU Make, tested with macOS and Linux. bfqutils is written in C and comes bundled with the Cloudflare fork of Zlib.

```sh
git clone https://github.com/<TODO>/bfqutils
cd bfqutils
make libz 
make release
make test     # optional
```

To dynamically link to a system Zlib instead of building a static library:

```sh
make z_dyn=1 release
```

## Citation

bfqutils hasn't appeared in any publication, but please do cite this repository if you find this software helpful in your research.

## bfqmerge

This tool is meant to be used as a drop-in replacement for [`fastp --merge`](https://github.com/OpenGene/fastp?tab=readme-ov-file#merge-paired-end-reads). The only differences from fastp are a smaller default minimum merge length and more aggressive base and quality score correction. (Additionally, bfqmerge seems to properly trim some polyG tails in cases where fastp fails to do so, though I don't quite understand why. Nevertheless, the polyG trimming algorithm isn't perfect and will occasionally fail.) Apart from that, bfqmerge is intended to be a lightweight replacement with extremely low memory usage (<2 MB) and very good single-thread performance (<15 min on 2x100M PE150 reads).

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
 -z        Compress the output as gzip.
 -q        Make the program quiet.
 -v        Print the version and exit.
 -h        Print this help message and exit.
 ```

### Implementation details

bfqmerge is a simple program, which performs the following operations:

1. Detect and trim polyG tails (`-g`). Failure to remove these can lead to bfqmerge thinking that sufficiently long stretches of Gs found in both reads of a pair are the overlapping parts of the reads.

2. Find the best overlap, depending on user settings (`-o`, `-d`, `-p`).

3. If no overlap is found, discard the reads. Otherwise, create a new merged read. For each overlapping position, the base and quality score is taken from whichever of the two reads has a better score in that position.

4. Check if the merged read passes quality filters (`-Q`, `-u`, `-n`). If yes, then compress the output if desired (`-z`), then write to `stdout`.

With the exception of `-o`, all default values are identical to fastp. These generally work quite well, though be aware that it is impossible to get read merging right 100% of the time. False positive merging events and false negative read discards will almost always occur given a sufficiently diverse set of reads.

