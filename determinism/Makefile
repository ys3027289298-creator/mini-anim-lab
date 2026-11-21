CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: determinism

determinism: src/determinism.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: determinism
	./determinism

clean:
	rm -f determinism
