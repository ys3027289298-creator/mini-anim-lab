#!/usr/bin/env bash
# End-to-end acceptance reproduction for the LOD/frustum coupling.
# Produces every deliverable under <outdir> (default: ./acceptance_out):
#   linked/   culled frames, contact sheet, triptychs, CSVs (linked pipeline)
#   legacy/   same, linkage off (the bug-reproducing pipeline)
#   hashes/   sha256 of two independent legacy runs, diffed for pixel-identity
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT="${1:-$ROOT/acceptance_out}"
FRAMES="${2:-240}"

echo "== build =="
make -C "$ROOT" clean all

echo "== self-test =="
"$ROOT/lod"

rm -rf "$OUT"
mkdir -p "$OUT/linked" "$OUT/legacy" "$OUT/legacy_b" "$OUT/hashes"

echo "== simulate linked=1 (coupled) =="
"$ROOT/lod" sim "$OUT/linked" 1 "$FRAMES"
echo "== simulate linked=0 (legacy), twice =="
"$ROOT/lod" sim "$OUT/legacy" 0 "$FRAMES"
"$ROOT/lod" sim "$OUT/legacy_b" 0 "$FRAMES"

echo "== pixel-identity check for linkage-off =="
(cd "$OUT/legacy" && find . -type f | sort | xargs sha256sum) > "$OUT/hashes/run_a.sha256"
(cd "$OUT/legacy_b" && find . -type f | sort | xargs sha256sum) > "$OUT/hashes/run_b.sha256"
diff "$OUT/hashes/run_a.sha256" "$OUT/hashes/run_b.sha256"
echo "legacy runs byte-identical: OK"
rm -rf "$OUT/legacy_b"

echo "== no-hole assertion =="
grep -q '^worst_linked_missing_px=0$' "$OUT/linked/run_meta.txt"
grep -q '^worst_linked_missing_px=0$' "$OUT/legacy/run_meta.txt"
echo "linked pipeline missing pixels vs LOD0 truth: 0"

echo "== PPM -> PNG =="
python3 -I "$HERE/ppm2png.py" "$OUT/linked"
python3 -I "$HERE/ppm2png.py" "$OUT/legacy"

echo
echo "deliverables in $OUT:"
echo "  frames/sheet/triptychs : PNG + PPM"
echo "  error_table.csv        : end-effector error, legacy vs welded"
echo "  switch_log.csv         : every LOD switch moment and its hole count"
echo "  frame_metrics.csv      : per-frame culling and coverage"
