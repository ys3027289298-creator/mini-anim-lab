CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: retarget

retarget: src/retarget.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: retarget
	./retarget

clean:
	rm -f retarget
