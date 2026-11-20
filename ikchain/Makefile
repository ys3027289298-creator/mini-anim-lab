CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ikchain

ikchain: src/ikchain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ikchain
	./ikchain

clean:
	rm -f ikchain
