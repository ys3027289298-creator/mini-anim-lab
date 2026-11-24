CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
PYTHON ?= python3 -I

all: momentum

momentum: src/momentum.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: momentum
	./momentum

# Regenerate the ten-trial CSVs under ../artifacts/momentum.
accept: momentum
	./momentum --accept ../artifacts/momentum

# Independent Python re-derivation, charts and REPORT.md from those CSVs.
verify:
	$(PYTHON) tools/verify_model.py ../artifacts/momentum
	$(PYTHON) tools/make_figures.py ../artifacts/momentum
	$(PYTHON) tools/report.py ../artifacts/momentum

clean:
	rm -f momentum
