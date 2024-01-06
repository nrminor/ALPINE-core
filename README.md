# ALPINE: Anachronistic Lineage and Persistent INfection Explorer...but in Rust
[![Open Source Starter Files](https://github.com/nrminor/ALPINE.rs/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/ALPINE.rs/actions/workflows/open-source-starter.yml) [![Rust CI](https://github.com/nrminor/ALPINE.rs/actions/workflows/rust-ci.yml/badge.svg)](https://github.com/nrminor/ALPINE.rs/actions/workflows/rust-ci.yml)

See [the original ALPINE workflow repo here](https://github.com/nrminor/ALPINE). This repo is currently little more than a personal project for learning Rust. However, it's possible that a Rust binary may eventually replace the core components of the ALPINE pipeline that are written in Julia. And in the long-term, this tool may even allow researchers to perform ALPINE-type analyses without invoking the entire Nextflow pipeline and associated Docker containers.

For now, though, certain ALPINE components _could_ be replaced with commands like the following:
```
# replace "-" symbols with "N" characters
alpine replace-gaps <FASTA>

alpine separate-by-month <FASTA> <METADATA>
```

To set up and try out the toolset, either install the Rust toolchain and build it yourself or run `./easy_install` in this directory after cloning it to your machine.
