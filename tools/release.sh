#!/usr/bin/env bash
# Raise the version to X.Y.Z in one commit: the workspace's in Cargo.toml, the
# exact aci-tallies-geom pin in crates/tallies/Cargo.toml (the Python package
# takes the workspace version) and Cargo.lock.
#
#   tools/release.sh 0.2.0
#
# Pushing that commit to master releases it: the Release workflow sees a
# version with no tag yet, tests, builds, publishes and tags vX.Y.Z.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

version="${1:-}"
[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.-]+)?$ ]] ||
  { echo "usage: tools/release.sh X.Y.Z" >&2; exit 2; }

[ -z "$(git status --porcelain)" ] ||
  { echo "the working tree has changes: commit or stash them first" >&2; exit 1; }
branch=$(git rev-parse --abbrev-ref HEAD)
[ "$branch" = master ] || echo "note: on branch $branch; releases are made from master" >&2

git fetch --quiet --tags origin 2>/dev/null || true
! git rev-parse -q --verify "refs/tags/v$version" >/dev/null ||
  { echo "v$version is already released: versions are never published twice" >&2; exit 1; }

old=$(sed -n '/^\[workspace.package\]/,/^\[/s/^version = "\(.*\)"/\1/p' Cargo.toml)
[ -n "$old" ] || { echo "no version in [workspace.package] of Cargo.toml" >&2; exit 1; }

# `sed -i.bak` works with both GNU and BSD sed.
sed -i.bak "/^\[workspace.package\]/,/^\[/s/^version = \".*\"/version = \"$version\"/" Cargo.toml
sed -i.bak "s/\(aci-tallies-geom = .*version = \"=\)[^\"]*\"/\1$version\"/" crates/tallies/Cargo.toml
rm -f Cargo.toml.bak crates/tallies/Cargo.toml.bak
grep -q "^version = \"$version\"" Cargo.toml ||
  { echo "failed to set the workspace version" >&2; exit 1; }
grep -q "aci-tallies-geom = .*version = \"=$version\"" crates/tallies/Cargo.toml ||
  { echo "failed to set the aci-tallies-geom pin" >&2; exit 1; }

# Only the workspace's own entries in Cargo.lock change.
cargo update --workspace --quiet

git commit --quiet -am "Release $version"
echo "Version $old -> $version committed in $(git rev-parse --short HEAD)."
echo "Pushing it to master releases it."
