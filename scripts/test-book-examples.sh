#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

examples_dir="docs/examples"
output_dir="target/book-examples"

cargo build
mkdir -p "$output_dir"

# Rust 1.100 moved the compiled artifacts to target/<profile>/build/<crate>/<hash>/out;
# older toolchains keep them in target/<profile>/deps. Search both. Paths must be
# absolute: rustc resolves relative --extern/-L paths against the source directory.
library_paths=(-L "dependency=$repo_root/target/debug/deps")
for dir in "$repo_root"/target/debug/build/*/*/out; do
    [ -d "$dir" ] || continue
    library_paths+=(-L "dependency=$dir")
done

for src in "$examples_dir"/*.rs; do
    name="$(basename "$src" .rs)"
    rustc \
        --edition=2021 \
        --extern "quanticsgrids=$repo_root/target/debug/libquanticsgrids.rlib" \
        "${library_paths[@]}" \
        "$src" \
        -o "$output_dir/$name"
    "$output_dir/$name"
done
