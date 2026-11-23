CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: curve

curve: src/curve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

ARTIFACT_DIR ?= $(CURDIR)/../artifacts/run_curve

test: curve
	mkdir -p $(ARTIFACT_DIR)
	ARTIFACT_DIR=$(ARTIFACT_DIR) ./curve

clean:
	rm -f curve
