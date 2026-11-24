#!/usr/bin/env bash
# Verify the unified pose format across skeleton / ik / curve / path.
#
# For each module (four independent passes):
#   1. build it with its own Makefile
#   2. new format: export -> import -> compare against the pre-export struct
#      (--pose-roundtrip), must print roundtrip=1
#   3. new format sample: a fresh --pose-export must match samples/<m>.pose
#   4. legacy format: --legacy-check against samples/legacy/<m>.csv must print
#      legacy=1 (old files still read back to the same internal state)
#   5. old default path: no-arg run must still match tests/samples/baseline/<m>.out
#
# Exit 0 only if every check for every module passes.

set -u

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 99
cd "$root" || exit 99

modules=(skeleton ik curve path)
tmp="$(mktemp -d)" || exit 99
trap 'rm -rf "$tmp"' EXIT

fail=0
pass=0
report() { # report <ok> <label>
    if [ "$1" -eq 1 ]; then
        printf '  PASS  %s\n' "$2"
        pass=$((pass + 1))
    else
        printf '  FAIL  %s\n' "$2"
        fail=$((fail + 1))
    fi
}

for m in "${modules[@]}"; do
    printf '\n== %s ==\n' "$m"

    if ! make -C "$m" >"$tmp/build.$m.log" 2>&1; then
        printf '  FAIL  build\n'
        sed 's/^/        /' "$tmp/build.$m.log"
        fail=$((fail + 1))
        continue
    fi
    report 1 "build"

    bin="./$m/$m"

    # 1. new format roundtrip against the in-memory struct
    out="$("$bin" --pose-roundtrip "$tmp/$m.roundtrip" 2>&1)"
    [ "$out" = "roundtrip=1" ] && report 1 "new-format roundtrip" \
                               || { report 0 "new-format roundtrip ($out)"; }

    # 2. exported bytes are stable and match the checked-in sample
    "$bin" --pose-export "$tmp/$m.pose" || true
    if diff -u "samples/$m.pose" "$tmp/$m.pose" >"$tmp/d.$m.pose"; then
        report 1 "new-format sample matches samples/$m.pose"
    else
        report 0 "new-format sample differs from samples/$m.pose"
        sed 's/^/        /' "$tmp/d.$m.pose"
    fi

    # 3. legacy file reads back to the same internal state
    out="$("$bin" --legacy-check "samples/legacy/$m.csv" 2>&1)"
    [ "$out" = "legacy=1" ] && report 1 "legacy read (samples/legacy/$m.csv)" \
                            || { report 0 "legacy read ($out)"; }

    # 4. the historical default output is byte-for-byte unchanged
    "$bin" >"$tmp/$m.default" 2>"$tmp/$m.default.err"
    if diff -u "tests/samples/baseline/$m.out" "$tmp/$m.default" >"$tmp/d.$m.default"; then
        report 1 "old default output unchanged"
    else
        report 0 "old default output changed"
        sed 's/^/        /' "$tmp/d.$m.default"
    fi
done

printf '\n==== summary ====\n'
printf 'pass=%d fail=%d\n' "$pass" "$fail"
if [ "$fail" -ne 0 ]; then
    printf 'FORMAT CHECK: FAILED\n'
    exit 1
fi
printf 'FORMAT CHECK: OK\n'
exit 0
