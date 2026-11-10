<!-- unified contribution policy: ../../docs/contributing/rust-repositories.md in tensor4all/tensor4all-meta -->
# Contributing to quanticsgrids-rs

The rules that apply to every Rust repository in `tensor4all` are in one place:
**[tensor4all-meta · docs/contributing/rust-repositories.md](https://github.com/tensor4all/tensor4all-meta/blob/main/docs/contributing/rust-repositories.md)**.
Read that first — it covers the issue-first workflow, what a pull request body has to contain, red-first tests, the
gates every repository runs, the review requirement, and the version policy.

This file only records what is specific to `quanticsgrids-rs`:

## Local facts

- **MSRV**: none is declared — no `rust-version` in `Cargo.toml` and no `rust-toolchain.toml`. CI pins `stable`, so
  there is no MSRV claim for CI to hold and no `msrv` job.
- **Feature combinations CI builds**: one. The single package declares no features, and `lint`, `test` and `docs` all
  run the default set.
- **Platforms claimed**: `ubuntu-22.04` for all three jobs. The release tests run through
  `cargo nextest run --release`, and the `docs` job additionally builds the mdBook and checks its executable examples
  with `./scripts/test-book-examples.sh`.
- **Local gate**: `cargo test --release`.

## Anything here that departs from the shared rules

One does: the `msrv` gate is not run, and no `rust-version` is declared for it to hold CI to. Declaring one and adding
the job is the way to close that gap; this file records the state and does not exempt the repository from the other
gates.
