CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: skeleton

skeleton: src/skeleton.c
	$(CC) $(CFLAGS) $< -o $@ -lm

ARTIFACT_DIR ?= $(CURDIR)/../artifacts/run_skeleton

test: skeleton
	mkdir -p $(ARTIFACT_DIR)
	ARTIFACT_DIR=$(ARTIFACT_DIR) ./skeleton

clean:
	rm -f skeleton
