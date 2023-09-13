# ALPINE: Anachronistic Lineage and Persistent INfection Explorer...but in Rust
See [the original ALPINE workflow repo here](https://github.com/nrminor/ALPINE). This is a personal project for Rust-learning purposes. It is unlikely that the workflow's components, which are written in Julia, R, and Python, will be replaced with the contents of this repo.

Still, the ALPINE components _could_ be replaced with commands like the following:
```
# replace "-" symbols with "N" characters
alpine replace-gaps <FASTA>

alpine separate-by-month <FASTA> <METADATA>
```
