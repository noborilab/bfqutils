# fqmerge: FastQ merging tool for PE reads

This tool is meant to be used as a drop-in replacement for `fastp --merge`. The only differences from fastp are a smaller default minimum merge length and more aggressive base and quality score correction. (Additionally, fqmerge seems to properly trim some polyG tails in cases where fastp fails to do so, though I don't quite understand why.) Apart from that, fqmerge is intended to be a lightweight replacement with extremely low memory usage and very performant single-thread usage.


