# Rust User Guide Design

## Goal

Add a Rust user guide that covers the same topics as the Julia documentation for
`QuanticsGrids.jl`, with executable examples and GitHub Pages publishing.

## Constraints

- The Rust crate is not yet published on crates.io.
- The guide must cover the same content as `../QuanticsGrids.jl/docs/src/index.md`.
- Example code should be executable and checked automatically.
- The guide should be published via GitHub Pages.

## Decisions

### Documentation Technology

Use `mdBook` for the user guide.

Rationale:

- It supports multi-page narrative documentation cleanly.
- It fits GitHub Pages deployment well.
- `mdbook test` can validate Rust code examples.
- It is better suited than `rustdoc` alone for a Julia-style user guide.

### Coverage

The Rust guide will cover the same content areas as the Julia guide:

1. Installation
2. Definition
3. Usage overview
4. Discretized grid
5. Creating a one-dimensional grid
6. Creating a d-dimensional grid
7. Inherent discrete grid
8. Using a base other than 2
9. Creating a function that takes a quantics index as input

The wording and examples will be adapted to the Rust API and idioms rather than
ported line-for-line from Julia syntax.

### Executable Examples

Primary examples in the guide will be written so they can be checked with
`mdbook test`.

Guidelines:

- Each code block should be self-contained.
- Each example should import the required crate items explicitly.
- Examples returning `Result` should use `fn main() -> quanticsgrids::Result<()>`.
- Floating-point examples should use tolerance-based assertions.

`cargo test --doc` will remain in place for API docs, but the user guide will be
validated separately through `mdbook test`.

### Repository Layout

Add:

- `book.toml`
- `docs/src/SUMMARY.md`
- `docs/src/index.md`
- `docs/src/installation.md`
- `docs/src/definition.md`
- `docs/src/discretized-grid.md`
- `docs/src/inherent-discrete-grid.md`
- `docs/src/non-binary-base.md`
- `docs/src/quantics-function.md`

Adjust:

- `README.md` to point to the GitHub Pages guide
- CI workflows to build and test the guide

### Installation Instructions

Because the crate is not yet published, installation examples will use the Git
dependency form:

```toml
[dependencies]
quanticsgrids = { git = "https://github.com/tensor4all/quanticsgrids-rs" }
```

Optionally mention:

```bash
cargo add quanticsgrids --git https://github.com/tensor4all/quanticsgrids-rs
```

### Publishing

Publish the generated `mdBook` site to GitHub Pages through a dedicated GitHub
Actions workflow.

CI should separately verify:

- `mdbook build`
- `mdbook test`

## Open Notes

- `Cargo.lock` should be removed from version control because this repository is
  a library crate.
- Existing license updates should be retained and finished alongside the docs
  work.
