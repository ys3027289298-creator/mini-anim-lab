CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: grapple

grapple: src/grapple.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: grapple
	./grapple

clean:
	rm -f grapple
