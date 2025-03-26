# fqmerge: FastQ merging tool for PE reads

This tool is meant to be used as a drop-in replacement for `fastp --merge`. The only differences from fastp are a smaller default minimum merge length and more aggressive base and quality score correction. (Additionally, fqmerge seems to properly trim some polyG tails in cases where fastp fails to do so, though I don't quite understand why.) Apart from that, fqmerge is intended to be a lightweight replacement with extremely low memory usage (<2 MB) and very performant single-thread runtime (<15 min on 100M read pairs).

## Installation

Requires gcc/clang and GNU Make, tested with macOS and Linux. After downloading:

```sh
make libz 
make release
make test     # optional
```

To dynamically link to Zlib instead of building a static library:

```sh
make z_dyn=1 release
```

## Usage

```
fqmerge v1.0  Copyright (C) 2025  Benjamin Jean-Marie Tremblay

Usage:  fqmerge [options] R1.fq[.gz] R2.fq[.gz] > merged.fq
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

