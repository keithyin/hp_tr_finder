# hp_tr_finder

[![CI](https://github.com/keithyin/hp_tr_finder/actions/workflows/ci.yml/badge.svg)](https://github.com/keithyin/hp_tr_finder/actions/workflows/ci.yml)
[![Release](https://github.com/keithyin/hp_tr_finder/actions/workflows/release.yml/badge.svg)](https://github.com/keithyin/hp_tr_finder/actions/workflows/release.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

homopolymer and tandem-repeats finder

## Install

```sh
cargo install hp_tr_finder --features cli
```

## Library

CLI dependencies (`clap`, `noodles-fasta`) are optional, behind the `cli` feature.
Depending on this crate as a library pulls in only `regex` and `intervaltree`:

```toml
[dependencies]
hp_tr_finder = "0.6.2"   # default features contain no CLI deps
```

## Build from source

```sh
cargo build --release --features cli
```
