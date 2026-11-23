CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lod

lod: src/lod.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lod
	./lod

debug:
	$(CC) -O0 -g -std=c11 -Wall -Wextra -DLOD_NO_MAIN tests/debug_frame.c -o /tmp/debug_frame -lm

clean:
	rm -f lod
