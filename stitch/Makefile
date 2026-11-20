CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: stitch

stitch: src/stitch.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: stitch
	./stitch

clean:
	rm -f stitch
