CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: coupling

coupling: src/coupling.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: coupling
	./coupling

clean:
	rm -f coupling
