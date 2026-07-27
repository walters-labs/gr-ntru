# GR-NTRU

Rust implementation of group ring NTRU over finite group rings.

Classical NTRU works in `Z[x]/(x^N - 1)`, which is the group ring `Z[C_N]`.
This crate replaces `C_N` by explicit cyclic, dihedral, and symmetric groups.
For noncommutative groups it uses the left-sided convention

```text
h = f^{-1} g mod q
e = p h r + m mod q
```

and decrypts by multiplying by `f` on the left.

The crate includes pure Rust group ring arithmetic and uses the
`fft-dihedral` and `fft-symmetric` crates, behind the default `fft` feature, for
fast multiplication and blockwise inversion when the field and group are compatible.

The Rust package and command are named `gr-ntru`; the library is imported as
`gr_ntru`.

## Install

```bash
cargo install gr-ntru
```

## Run

```bash
cargo run -- --profile quick
cargo run -- --profile fft

# After installation:
gr-ntru --profile quick
```

## Test

```bash
cargo test
cargo test --no-default-features
```

## Status

This is research software, not production cryptography.
