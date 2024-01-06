just default:
    just --list

dev-rs:
    cargo install scidataflow
    cargo install nu --features dataframe

macos:
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    brew install pre-commit

alpine-rs:
    cargo build --release
