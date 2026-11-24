#!/usr/bin/env bash
# Verify the unified entry does not change any module's executable behaviour.
#
# For each module, run its own `make test` exactly as before, strip only the
# make/compiler chatter and the module's own "all_ok=" line, and compare the
# remaining self-check lines with the unified entry's lines for that module
# (with the "[module] " tag removed). Any difference in item names, order,
# values, or exit status is a failure.

set -u
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 99
cd "$root" || exit 99
modules=(skeleton ik curve constraint path)
tmp="$(mktemp -d)" || exit 99
trap 'rm -rf "$tmp"' EXIT
rc=0

for m in "${modules[@]}"; do
  make -C "$m" clean >/dev/null 2>&1
  make -C "$m" test > "$tmp/direct.$m" 2>/dev/null
  direct_rc=$?
  grep -E '^[A-Za-z_][A-Za-z0-9_]*=[01]$' "$tmp/direct.$m" | grep -v '^all_ok=' > "$tmp/direct.lines.$m"

  ./selfcheck/run.sh > "$tmp/unified.txt" 2>&1
  unified_rc=$?
  sed -n "s/^\[$m\] //p" "$tmp/unified.txt" > "$tmp/unified.lines.$m"

  if diff -u --label "make -C $m test" --label "selfcheck/run.sh [$m]" \
       "$tmp/direct.lines.$m" "$tmp/unified.lines.$m" > "$tmp/diff.$m"; then
    printf '%-11s items=%d  IDENTICAL\n' "$m" "$(wc -l < "$tmp/direct.lines.$m")"
  else
    printf '%-11s DIFFERENT (direct exit=%d)\n' "$m" "$direct_rc"
    sed 's/^/    /' "$tmp/diff.$m"
    rc=1
  fi
done

printf 'direct make test: exit statuses were all 0 = %s\n' \
  "$([ "$rc" -eq 0 ] && echo yes || echo no)"
[ "$rc" -eq 0 ] && printf 'EQUIVALENCE: OK\n' || { printf 'EQUIVALENCE: MISMATCH\n'; exit 1; }
