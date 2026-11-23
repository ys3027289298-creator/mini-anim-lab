#!/usr/bin/env bash
# Full acceptance run. Every check required by the task:
#   1. four modules, standalone self-test, unified + legacy paths
#   2. legacy stdout byte-identical to the pre-change baseline (entry off)
#   3. illegal-parameter cases rejected with the owning module named
#   4. cross-module pose consistency through the unified entry
#   5. per-module before/after pose numbers and two-tier timing
# Outputs are archived under acceptance/results/ and summarized in
# acceptance/REPORT.md.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RES="$ROOT/acceptance/results"
mkdir -p "$RES"
MODULES="skeleton ik curve constraint"
overall=0

run() { echo "+ $*"; "$@"; }

echo "=== 1/5 standalone module self-tests (both paths) ==="
for m in $MODULES; do
  log="$RES/${m}_selftest.log"
  if "$ROOT/acceptance/verify_module.sh" "$m" > "$log" 2>&1; then
    echo "$m: PASS (log: acceptance/results/${m}_selftest.log)"
  else
    echo "$m: FAIL"; cat "$log"; overall=1
  fi
done

echo "=== 2/5 baseline byte-identity (unified entry OFF) ==="
ident=1
for m in $MODULES; do
  if diff -q "$ROOT/acceptance/baseline_out/$m.stdout.txt" \
             "$ROOT/acceptance/legacy_out/$m.txt" >/dev/null; then
    echo "$m: byte-identical to baseline"
  else
    echo "$m: DIFFERS"; ident=0; overall=1
  fi
done
echo "$ident" > "$RES/baseline_identity.status"

echo "=== 3/5 illegal parameters name the module ==="
if make -C "$ROOT/acceptance" -s clean test > "$RES/invalid.log" 2>&1; then
  echo "invalid: PASS ($(grep -c ' ok=1' "$RES/invalid.log") cases)"
else
  echo "invalid: FAIL"; cat "$RES/invalid.log"; overall=1
fi

echo "=== 4/5 cross-module pose consistency ==="
if make -C "$ROOT/pose" -s clean test > "$RES/pose.log" 2>&1; then
  echo "pose: PASS"; grep -E 'pose_|poses_consistent' "$RES/pose.log"
else
  echo "pose: FAIL"; cat "$RES/pose.log"; overall=1
fi

echo "=== 5/5 pose snapshots + two-tier timing + before/after values ==="
if "$ROOT/acceptance/compare_poses.sh" > "$RES/compare.log" 2>&1; then
  echo "pose before/after: PASS"; cat "$RES/compare.log"
else
  echo "pose before/after: FAIL"; cat "$RES/compare.log"; overall=1
fi
if make -C "$ROOT/acceptance/bench" -s clean test > "$RES/bench.log" 2>&1; then
  echo "bench: PASS"; grep -E 'time_|byte_identical|match=0' "$RES/bench.log"
  if grep -q 'match=0' "$RES/bench.log" || grep -q 'poses_byte_identical=0' "$RES/bench.log"; then
    echo "pose mismatch found in bench"; overall=1
  fi
else
  echo "bench: FAIL"; cat "$RES/bench.log"; overall=1
fi

echo
if [[ $overall -eq 0 ]]; then
  echo "ALL ACCEPTANCE CHECKS PASSED"
else
  echo "ACCEPTANCE FAILURES PRESENT"
fi
exit $overall
