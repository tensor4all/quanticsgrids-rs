#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

examples_dir="docs/examples"
output_dir="target/book-examples"

cargo build
mkdir -p "$output_dir"

for src in "$examples_dir"/*.rs; do
    name="$(basename "$src" .rs)"
    rustc \
        --edition=2021 \
        --extern "quanticsgrids=target/debug/libquanticsgrids.rlib" \
        -L dependency=target/debug/deps \
        "$src" \
        -o "$output_dir/$name"
    "$output_dir/$name"
done
