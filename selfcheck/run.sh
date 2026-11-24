#!/usr/bin/env bash
# Unified self-check entry for mini-anim-lab.
#
# Runs each module's own build + self-check binary and merges their results.
# It does not modify module sources, Makefiles, or executable behaviour: each
# module is built with its own Makefile and its binary is run as-is. Only its
# stdout (the "name=1"/"name=0" lines) is parsed.
#
# Every result line is tagged with its origin module. Each module's expected
# check names are read from its source's check("...") calls, so an item that
# fails to run (crash, early return) is reported as MISSING rather than
# silently dropped.
#
# Exit status: 0 only if every expected item ran and passed.

set -u

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 99
cd "$root" || exit 99

modules=(skeleton ik curve constraint path)
tmp="$(mktemp -d)" || exit 99
trap 'rm -rf "$tmp"' EXIT

total_pass=0
total_fail=0
total_missing=0
failed_modules=0
summary=""

for m in "${modules[@]}"; do
  src="$m/src/$m.c"
  bin="$m/$m"
  printf '\n== %s ==\n' "$m"

  if [ ! -f "$src" ]; then
    printf '[%s] SOURCE MISSING: %s\n' "$m" "$src"
    failed_modules=$((failed_modules + 1))
    summary+="$m pass=0 fail=0 missing=0 (source missing)"$'\n'
    continue
  fi

  # Expected items: literal check("name", ...) calls in the module source.
  grep -oE 'check\("[A-Za-z_][A-Za-z0-9_]*"' "$src" \
    | sed -E 's/^check\("//; s/"$//' > "$tmp/expected.$m"

  # Build with the module's own Makefile (default target).
  if ! make -C "$m" > "$tmp/build.$m.log" 2>&1; then
    printf '[%s] BUILD FAILED\n' "$m"
    sed 's/^/    /' "$tmp/build.$m.log"
    failed_modules=$((failed_modules + 1))
    summary+="$m pass=0 fail=0 missing=0 (build failed)"$'\n'
    continue
  fi

  if [ ! -x "$bin" ]; then
    printf '[%s] BINARY MISSING: %s\n' "$m" "$bin"
    failed_modules=$((failed_modules + 1))
    summary+="$m pass=0 fail=0 missing=0 (binary missing)"$'\n'
    continue
  fi

  # Run the module's self-check binary exactly as `make test` would.
  "$bin" > "$tmp/out.$m" 2> "$tmp/err.$m"
  run_rc=$?

  pass=0
  fail=0
  missing=0
  all_ok="-"
  declare -A seen=()

  while IFS= read -r line; do
    if [[ $line =~ ^([A-Za-z_][A-Za-z0-9_]*)=([01])$ ]]; then
      name="${BASH_REMATCH[1]}"
      val="${BASH_REMATCH[2]}"
      if [ "$name" = "all_ok" ]; then
        all_ok="$val"
        continue
      fi
      seen["$name"]=1
      printf '[%s] %s=%s\n' "$m" "$name" "$val"
      if [ "$val" = "1" ]; then
        pass=$((pass + 1))
      else
        fail=$((fail + 1))
      fi
    else
      printf '[%s] (unparsed output) %s\n' "$m" "$line"
    fi
  done < "$tmp/out.$m"

  # Any item the source declares but that never reported is a lost item.
  while IFS= read -r name; do
    [ -n "$name" ] || continue
    if [ -z "${seen[$name]+x}" ]; then
      printf '[%s] MISSING %s\n' "$m" "$name"
      missing=$((missing + 1))
    fi
  done < "$tmp/expected.$m"

  unset seen

  # Guard rails for "no items lost" and for the modules' own summary line.
  expected_count=$(wc -l < "$tmp/expected.$m")
  if [ "$expected_count" -eq 0 ]; then
    printf '[%s] NO EXPECTED ITEMS FOUND IN SOURCE\n' "$m"
    missing=$((missing + 1))
  fi
  if [ "$run_rc" -ne 0 ] && [ "$fail" -eq 0 ] && [ "$missing" -eq 0 ]; then
    printf '[%s] ANOMALY: exited %d with no failing item reported\n' "$m" "$run_rc"
    fail=$((fail + 1))
  fi
  if [ "$fail" -eq 0 ] && [ "$missing" -eq 0 ]; then want_all_ok=1; else want_all_ok=0; fi
  if [ "$all_ok" != "-" ] && [ "$all_ok" != "$want_all_ok" ]; then
    printf '[%s] ANOMALY: all_ok=%s disagrees with item results\n' "$m" "$all_ok"
    fail=$((fail + 1))
  fi
  if [ -s "$tmp/err.$m" ]; then
    printf '[%s] stderr:\n' "$m"
    sed 's/^/    /' "$tmp/err.$m"
  fi

  printf -- '-- %s: pass=%d fail=%d missing=%d\n' "$m" "$pass" "$fail" "$missing"
  summary+="$m pass=$pass fail=$fail missing=$missing"$'\n'

  total_pass=$((total_pass + pass))
  total_fail=$((total_fail + fail))
  total_missing=$((total_missing + missing))
  if [ $((fail + missing)) -gt 0 ]; then
    failed_modules=$((failed_modules + 1))
  fi
done

printf '\n==== summary ====\n'
printf '%s' "$summary"
printf -- '-----------------\n'
printf 'total: pass=%d fail=%d missing=%d modules_failed=%d\n' \
  "$total_pass" "$total_fail" "$total_missing" "$failed_modules"

if [ $((total_fail + total_missing)) -gt 0 ] || [ "$failed_modules" -gt 0 ]; then
  printf 'SELFCHECK: FAILED\n'
  exit 1
fi
printf 'SELFCHECK: OK\n'
exit 0
