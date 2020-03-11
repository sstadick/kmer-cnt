# kmer-cnt

These are rust clones of https://github.com/lh3/kmer-cnt.

## Input data

```
wget https://github.com/lh3/kmer-cnt/releases/download/v0.1/M_abscessus_HiSeq_10M.fa.gz
```

## Build and run

In each directory, you can run `cargo build --release` and a binary will
be created in the target/release dir with the name indicated in the
Cargo.toml file.

All benchmarking is being done via gzcat M_abscessus_HiSeq_10M.fa.gz
into each of these test programs.

Fasta file parseing may be different than the original versions.
