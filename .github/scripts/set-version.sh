#!/usr/bin/env bash
# Stamp the release version from the pushed tag (vX.Y.Z) into the Cargo
# manifests, so the tag is the single source of truth for the version.
# pyproject.toml takes its version from Cargo.toml (dynamic = ["version"]).
# Uses `sed -i.bak` so it works with both GNU sed (Linux, Windows Git Bash)
# and BSD sed (macOS).
set -euo pipefail

version="${GITHUB_REF_NAME#v}"
if ! [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+([-+][0-9A-Za-z.-]+)?$ ]]; then
  echo "::error::tag '$GITHUB_REF_NAME' is not of the form vX.Y.Z"; exit 1
fi

sed -i.bak "/^\[workspace.package\]/,/^\[/s/^version = \".*\"/version = \"$version\"/" Cargo.toml
sed -i.bak "s/\(aci-tallies-geom = .*version = \"=\)[^\"]*\"/\1$version\"/" crates/tallies/Cargo.toml
rm -f Cargo.toml.bak crates/tallies/Cargo.toml.bak

grep -q "^version = \"$version\"" Cargo.toml \
  || { echo "::error::failed to set version in Cargo.toml"; exit 1; }
grep -q "aci-tallies-geom = .*version = \"=$version\"" crates/tallies/Cargo.toml \
  || { echo "::error::failed to set aci-tallies-geom pin in crates/tallies/Cargo.toml"; exit 1; }
echo "Release version: $version"
