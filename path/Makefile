CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: path

path: src/path.c
	$(CC) $(CFLAGS) $< -o $@ -lm

ARTIFACT_DIR ?= $(CURDIR)/../artifacts/run_path

test: path
	mkdir -p $(ARTIFACT_DIR)
	ARTIFACT_DIR=$(ARTIFACT_DIR) ./path

clean:
	rm -f path
