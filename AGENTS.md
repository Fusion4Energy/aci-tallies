# aci-tallies

Readers for MCNP tally output (meshtal, MCTAL, D1SUNED CuV) in Rust, with
Python bindings. `README.md` is the user-facing reference; this file is for
whoever changes the code.

## Priorities

In this order, for every change:

1. **Correctness.** A value is right or it is reported: never a guess. What
   can't be read becomes a `failed` section with its reason, or NaN with a
   `Warning` that has the line number. Unknown input (a column label, a
   keyword, a format) is refused by name, not silently mapped to something
   close. MCNP itself is the reference for what a file means: when in doubt,
   run MCNP on a tiny input and write down what it did (in a test or the
   README), rather than reading it off another tool.
2. **Performance.** Files can be larger than memory and machines have many
   cores. Every reader follows the same shape: `common::source::Source`, one
   parallel `scan::index` pass that marks the lines it needs, then reads only
   the bytes asked for, parsed in parallel with rayon. Memory scales with the
   result, not the file. Measure (time and peak RSS) before and after a change
   to a reader; `examples/` time a file.
3. **Simplicity.** One copy of each thing: vector maths in
   `aci_tallies_geom`, number parsing in `common/float.rs`, line scanning in
   `common/scan.rs`. A new reader reuses them rather than growing its own. No
   public API without a user for it; what a caller can derive from what is
   already exposed stays out.

## Working on it

- The checks CI runs (`.github/workflows/ci.yml`), all must pass:
  ```sh
  cargo fmt --all --check
  cargo clippy --workspace --all-targets -- -D warnings
  cargo test --workspace
  pytest tests/python && mypy python/aci_tallies tests/python && ruff check python tests
  ```
  Build the extension with `maturin develop --release` (debug is 10–30x
  slower). It writes `python/aci_tallies/_native*.so` into the tree, which
  every editable install of the repo then loads.
- A change to a reader is checked for regressions on real files, not only
  the fixtures: read every tally of every file before and after and compare
  values, `failed` lists and warnings. Real output beats synthetic fixtures:
  D1SUNED's own files have shown forms the fixtures didn't.
- Every reader's parse must not depend on how the file is cut into chunks:
  test with small chunk sizes (`build(source, ..., chunk)`).
- Fixtures under `tests/fixtures` are MCNP output; each `.mcnp` says how to
  regenerate it.
- Docs and comments go stale with the code: update the README and doc
  comments in the same change.

## Versions and releases

The version lives in `Cargo.toml` (`[workspace.package]`, and the exact
`aci-tallies-geom` pin in `crates/tallies/Cargo.toml`). A push to `master`
whose version has no tag yet is released
(`.github/workflows/release.yml`); nobody types a tag. So:

- Feature branches and PRs never change the version.
- To release, on a clean, up-to-date `master`: `tools/release.sh X.Y.Z`
  (sets both, updates `Cargo.lock`, commits), then push. The workflow tests,
  builds the wheels, publishes to PyPI and crates.io, and the GitHub release
  creates the tag `vX.Y.Z` on that commit.
- Versions are never published twice.
