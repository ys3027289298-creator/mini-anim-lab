CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sync

sync: src/sync.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sync
	./sync

clean:
	rm -f sync
