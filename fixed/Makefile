CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: fixed

fixed: src/fixed.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: fixed
	./fixed

clean:
	rm -f fixed
