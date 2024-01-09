# ALPINE (Anachronistic Lineage and Persistent INfection Explorer) Core Utilities
[![Open Source Starter Files](https://github.com/nrminor/ALPINE.rs/actions/workflows/open-source-starter.yml/badge.svg)](https://github.com/nrminor/ALPINE.rs/actions/workflows/open-source-starter.yml) [![Rust CI](https://github.com/nrminor/ALPINE.rs/actions/workflows/rust-ci.yml/badge.svg)](https://github.com/nrminor/ALPINE.rs/actions/workflows/rust-ci.yml)

See [the original ALPINE workflow repo here](https://github.com/nrminor/ALPINE). This repo contains the in-development Rust source code meant to replace the core components of the ALPINE pipeline that were previously written in Julia. This tool has a rich command line interface intended to allow researchers to perform ALPINE-type analyses without invoking the entire Nextflow pipeline and associated Docker containers.

When you run `alpine --help`, the tool usage that comes up looks like this:
```

_____/\\\\\\\\\_____/\\\______________/\\\\\\\\\\\\\____/\\\\\\\\\\\__/\\\\\_____/\\\__/\\\\\\\\\\\\\\\_
 ___/\\\\\\\\\\\\\__\/\\\_____________\/\\\/////////\\\_\/////\\\///__\/\\\\\\___\/\\\_\/\\\///////////__
  __/\\\/////////\\\_\/\\\_____________\/\\\_______\/\\\_____\/\\\_____\/\\\/\\\__\/\\\_\/\\\_____________
   _\/\\\_______\/\\\_\/\\\_____________\/\\\\\\\\\\\\\/______\/\\\_____\/\\\//\\\_\/\\\_\/\\\\\\\\\\\_____
    _\/\\\\\\\\\\\\\\\_\/\\\_____________\/\\\/////////________\/\\\_____\/\\\\//\\\\/\\\_\/\\\///////______
     _\/\\\/////////\\\_\/\\\_____________\/\\\_________________\/\\\_____\/\\\_\//\\\/\\\_\/\\\_____________
      _\/\\\_______\/\\\_\/\\\_____________\/\\\_________________\/\\\_____\/\\\__\//\\\\\\_\/\\\_____________
       _\/\\\_______\/\\\_\/\\\\\\\\\\\\\\\_\/\\\______________/\\\\\\\\\\\_\/\\\___\//\\\\\_\/\\\\\\\\\\\\\\\_
        _\///________\///__\///////////////__\///______________\///////////__\///_____\/////__\///////////////__

ALPINE: Anachronistic Lineage and Persistent INfection Explorer
===============================================================

Command line interface for the core Rust components of ALPINE.
These commands are called in the full pipeline, which is imple-
mented in Nextflow alongside `seqkit`, `csvtk`, `nushell`,
`vsearch`, and a bin of bespoke Python scripts. However, users
may also use the commands available in this crate to run
similar analyses themselves via the command line.


Usage: alpine [OPTIONS] [COMMAND]

Commands:
  replace-gaps       Replace FASTA gap symbols '-' with masked bases 'N'.
  filter-by-n        Filter out FASTA records with more than the desired number of masked 'N' bases.
  separate-by-month  Use collection dates from FASTA record metadata to sort all FASTA records into a separate FASTA for each year-month combination.
  distance-matrix    Compute a symmetric pairwise distance matrix based on how dissimilar sequences in the provided FASTA are to one another.
  help               Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose...         Increase logging verbosity
  -q, --quiet...           Decrease logging verbosity
  -t, --threads <THREADS>  [default: 3]
  -h, --help               Print help
```

To set up and try out the toolset, either install the Rust toolchain and build it yourself or run `./easy_install` in this directory after cloning it to your machine. Once the ALPINE core utilities are ready for the bigtime, they will be made available via [crates.io](https://crates.io/).
