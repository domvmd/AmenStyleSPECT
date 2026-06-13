#!/usr/bin/env bash
# Build a distributable zip of the AmenStyleSPECT extension.
# Usage:  ./AmenStyleSPECT/package.sh   (run from anywhere)
# Output: <project root>/AmenStyleSPECT-<version>.zip
set -euo pipefail

VERSION="0.6.0"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../AmenStyleSPECT
ROOT="$(dirname "$HERE")"                               # project root (parent)
NAME="AmenStyleSPECT"
OUT="$ROOT/${NAME}-${VERSION}.zip"

rm -f "$OUT"
cd "$ROOT"
# Exclude editor/OS junk and Python bytecode caches.
zip -r "$OUT" "$NAME" \
    -x "*.DS_Store" \
    -x "*__pycache__*" \
    -x "*.pyc" \
    -x "*/.git/*" >/dev/null

echo "Wrote: $OUT"
echo "Contents:"
unzip -l "$OUT"
