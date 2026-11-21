CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint unified

constraint: src/constraint.c
	$(CC) $(CFLAGS) $< -o $@ -lm

unified: src/unified.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: constraint unified
	./constraint
	./unified

clean:
	rm -f constraint unified
