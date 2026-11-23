CC ?= cc
PYTHON ?= python3

MODULES := skeleton ik curve constraint path

.PHONY: all test test-% report clean help

all: test

# One command, all five modules, all five evidence kinds. Non-zero exit if any
# module fails or any evidence is missing; artifacts/selfcheck.{json,txt}.
test:
	$(PYTHON) -I tools/selfcheck.py --cc $(CC)

# Same runner, one module at a time.
test-%:
	$(PYTHON) -I tools/selfcheck.py --cc $(CC) --module $*

# Re-print the last recorded result without rebuilding.
report:
	@test -f artifacts/selfcheck.txt && cat artifacts/selfcheck.txt || \
		{ echo "no report yet; run 'make test' first"; exit 1; }

# Each module's own self-check, unchanged.
per-module:
	@set -e; for m in $(MODULES); do $(MAKE) -C $$m test; done

clean:
	@for m in $(MODULES); do $(MAKE) -C $$m clean; done
	rm -rf artifacts/build artifacts/selfcheck.json artifacts/selfcheck.txt

help:
	@echo "make test        unified self-check for all five modules"
	@echo "make test-<mod>  unified self-check for one module"
	@echo "make per-module  original per-directory 'make test'"
	@echo "make report      print the last recorded evidence"
	@echo "make clean       remove build outputs and recorded evidence"
