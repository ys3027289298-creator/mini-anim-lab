#!/usr/bin/env bash
# Verify one module standalone, in both call paths:
#   unified (MAL_UNIFIED=1, default) and legacy (MAL_UNIFIED=0).
# Legacy stdout is diffed byte-for-byte against the archived baseline.
# Usage: verify_module.sh skeleton|ik|curve|constraint
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mod="${1:?usage: verify_module.sh <skeleton|ik|curve|constraint>}"
case "$mod" in
  skeleton|ik|curve|constraint) ;;
  *) echo "unknown module: $mod" >&2; exit 2;;
esac

fail=0
echo "### module: $mod"

echo "-- unified build (mal_solve entry) --"
make -C "$ROOT/$mod" -s clean test || fail=1

echo "-- legacy build (direct static routines) --"
make -C "$ROOT/$mod" -s clean all MAL_UNIFIED=0 >/dev/null
"$ROOT/$mod/$mod" > "$ROOT/acceptance/legacy_out/$mod.txt" 2>&1 || fail=1
cat "$ROOT/acceptance/legacy_out/$mod.txt"

echo "-- byte-identical check vs archived baseline --"
if diff -u "$ROOT/acceptance/baseline_out/$mod.stdout.txt" \
           "$ROOT/acceptance/legacy_out/$mod.txt"; then
  echo "baseline_match=IDENTICAL"
else
  echo "baseline_match=DIFFERS"; fail=1
fi

make -C "$ROOT/$mod" -s clean all >/dev/null   # restore default build
echo "result=$([[ $fail -eq 0 ]] && echo PASS || echo FAIL)"
exit $fail
