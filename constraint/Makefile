CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint

constraint: src/constraint.c
	$(CC) $(CFLAGS) $< -o $@ -lm

ARTIFACT_DIR ?= $(CURDIR)/../artifacts/run_constraint

test: constraint
	mkdir -p $(ARTIFACT_DIR)
	ARTIFACT_DIR=$(ARTIFACT_DIR) ./constraint

clean:
	rm -f constraint
