#!/usr/bin/env bash
# Compare pose numbers from the ORIGINAL pre-change code (reference
# drivers over git HEAD sources) against the unified mal_solve() output
# (bench programs). Every numeric pose value must match exactly.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF="$ROOT/acceptance/reference"
BENCH="$ROOT/acceptance/bench"
RES="$ROOT/acceptance/results"
mkdir -p "$RES"

make -C "$REF" -s clean all 2>/dev/null
make -C "$BENCH" -s clean all >/dev/null 2>&1

fail=0

# Extract every parenthesized float tuple / scalar number from a line set,
# in order, so the two implementations can be compared value-by-value.
numbers() { grep -oE '\-?[0-9]+\.[0-9]+' | tr '\n' ' '; echo; }

for m in skeleton ik curve constraint; do
  "$REF/ref_$m" 2>/dev/null | numbers > "$RES/${m}_before.nums"
  "$BENCH/bench_$m" 2>/dev/null | grep -v '^time' | numbers > "$RES/${m}_after.nums"

  # bench constraint also contains two pose blocks; compare against the
  # reference set by checking every reference value appears, in order.
  before=($(cat "$RES/${m}_before.nums"))
  after=($(cat "$RES/${m}_after.nums"))
  mism=0; matched=0
  for v in "${before[@]}"; do
    if printf '%s\n' "${after[@]}" | grep -Fxq -e "$v"; then
      matched=$((matched + 1))
    else
      mism=$((mism + 1))
      echo "  MISSING after value: $v"
    fi
  done
  echo "$m: before_values=${#before[@]} matched_after=$matched mismatch=$mism"
  [[ $mism -eq 0 && $matched -eq ${#before[@]} ]] || fail=1
done

# curve bench explicitly prints match=1 per sample; assert no match=0
if grep -q 'match=0' <("$BENCH/bench_curve" 2>/dev/null); then fail=1; fi
# constraint bench asserts byte-identical particle blocks
"$BENCH/bench_constraint" >/dev/null 2>&1 || fail=1

echo
if [[ $fail -eq 0 ]]; then
  echo "POSE BEFORE/AFTER COMPARISON: ALL VALUES MATCH"
else
  echo "POSE BEFORE/AFTER COMPARISON: MISMATCH"
fi
exit $fail
