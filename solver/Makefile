CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
PYTHON ?= python3 -I

all: solver

solver: src/solver.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: solver
	./solver

# Regenerate every acceptance artefact under ../artifacts/solver.
accept: solver
	./solver --accept ../artifacts/solver

# Independent cross-check: re-rasterise the pose CSV in Python, count red PPM
# pixels, render review PNG/GIF, then build REPORT.md. Fails if fixed-order
# penetration is not zero.
verify:
	$(PYTHON) tools/verify_pixels.py ../artifacts/solver
	$(PYTHON) tools/make_frames.py ../artifacts/solver
	$(PYTHON) tools/report.py ../artifacts/solver

clean:
	rm -f solver
