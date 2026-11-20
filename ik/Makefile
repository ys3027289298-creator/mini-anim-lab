CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ik

ik: src/ik.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ik
	./ik

clean:
	rm -f ik
