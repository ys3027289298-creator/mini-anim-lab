CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ik

ik: src/ik.c
	$(CC) $(CFLAGS) $< -o $@ -lm

ARTIFACT_DIR ?= $(CURDIR)/../artifacts/run_ik

test: ik
	mkdir -p $(ARTIFACT_DIR)
	ARTIFACT_DIR=$(ARTIFACT_DIR) ./ik

clean:
	rm -f ik
