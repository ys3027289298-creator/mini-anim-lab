#!/usr/bin/env bash
# End-to-end reproduction for the momentum module.
#
#   ./reproduce.sh [OUTDIR]
#
# Builds, runs the C self-test, regenerates the ten-trial CSVs, re-derives every
# number in Python and rewrites the figures and REPORT.md. Exits non-zero if the
# C self-test or the independent Python check fails.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE="$(cd "$HERE/.." && pwd)"
ROOT="$(cd "$MODULE/.." && pwd)"
OUT="${1:-$ROOT/artifacts/momentum}"

echo "== build =="
make -C "$MODULE" clean all

echo "== C self-test =="
"$MODULE/momentum"

echo "== generate CSVs -> $OUT =="
mkdir -p "$OUT"
"$MODULE/momentum" --accept "$OUT"

echo "== independent Python re-derivation =="
python3 -I "$HERE/verify_model.py" "$OUT"

echo "== figures + report =="
python3 "$HERE/make_figures.py" "$OUT"
python3 -I "$HERE/report.py" "$OUT"

echo "== done: $OUT =="
