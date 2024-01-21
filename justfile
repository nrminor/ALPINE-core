just default:
    just --list

macos:
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install pre-commit

dev-rs:
    cargo install scidataflow
    cargo install nu --features dataframe

alpine-core:
    cargo build --path .
    cargo clean
alias ac := alpine-core

all:
    just macos
    just dev-rs
    just ac
