CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: path

path: src/path.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: path
	./path

clean:
	rm -f path
