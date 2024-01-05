# ALPINE: Anachronistic Lineage and Persistent INfection Explorer...but in Rust
See [the original ALPINE workflow repo here](https://github.com/nrminor/ALPINE). This repo is currently little more than a personal project for learning Rust. However, it's possible that a Rust binary may eventually replace the core components of the ALPINE pipeline that are written in Julia. And in the long-term, this tool may even allow researchers to perform ALPINE-type analyses without invoking the entire Nextflow pipeline and associated Docker containers.

For now, though, certain ALPINE components _could_ be replaced with commands like the following:
```
# replace "-" symbols with "N" characters
alpine replace-gaps <FASTA>

alpine separate-by-month <FASTA> <METADATA>
```
