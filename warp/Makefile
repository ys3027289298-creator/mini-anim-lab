CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: warp

warp: src/warp.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: warp
	./warp

clean:
	rm -f warp
