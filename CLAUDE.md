# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**bfqutils** is a C toolkit for FastQ sequencing data processing. It bundles three subcommands into one binary:
- `merge` — paired-end read merging (low-memory streaming, <2 MB)
- `trimse` — single-end adapter trimming and quality filtering
- `stats` — FastQ quality control and statistics

## Build

```bash
make libz        # compile bundled Cloudflare zlib fork as static lib (first time only)
make             # release build with -O3 (default target: release)
make debug       # debug build with -g3, AddressSanitizer + UBSan enabled
make z_dyn=1     # link against system zlib instead of bundled static lib
make clean       # full clean (libz + object files + binary)
```

The build requires only gcc/clang and GNU Make.

## Tests

```bash
make test        # builds release then runs test/test.sh
```

`test/test.sh` diffs the output of `bfqutils merge` and `bfqutils trimse` against reference files in `test/`. There is no per-tool test isolation; running `make test` exercises both tools end-to-end.

Reference files: `test/merged.fq`, `test/testR1.trim.fq`.

## Architecture

```
src/bfqutils.c       — main(); dispatches subcommand to main_merge/main_trimse/main_stats
src/bfqmerge.c       — merge implementation
src/bfqtrimse.c      — single-end trim implementation
src/bfqstats.c       — stats implementation
src/bfqmerge.h       — declares int main_merge(int, char **)
src/bfqtrimse.h      — declares int main_trimse(int, char **)
src/bfqstats.h       — declares int main_stats(int, char **)
src/version.h        — BFQUTILS_VERSION / BFQUTILS_YEAR constants
src/kseq.h           — vendored MIT FASTA/FASTQ streaming reader (Attractive Chaos)
libs/zlib/           — vendored Cloudflare zlib fork; compiled to libs/zlib/libz.a
```

Each tool is entirely self-contained in its `.c` file. There are no inter-tool dependencies. All three share the same streaming pattern via `kseq.h` and zlib for gzip I/O.

## Key implementation details

- **Standard:** GNU C99 (`-std=gnu99`). No C11/C17 features.
- **Quality encoding:** PHRED+33 assumed throughout; no auto-detection.
- **stdin support:** All tools accept `-` as the input filename for stdin.
- **Output:** Always written to stdout; optional gzip compression via `-z`.
- **Zlib:** Compiled as `libs/zlib/libz.a` and statically linked by default. The `configure.log` and compiled artifacts in `libs/zlib/` are untracked — run `make libz` to (re)generate them.
- **No threads:** Deliberately single-threaded; performance target is ~10–15 min for 200 M PE reads.
- **Adding a new subcommand:** create `src/bfqNEW.c` + `src/bfqNEW.h`, add its object to the Makefile, and add a dispatch branch in `src/bfqutils.c`.
