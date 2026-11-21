CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint

constraint: src/constraint.c src/chain.c src/chain.h
	$(CC) $(CFLAGS) src/constraint.c src/chain.c -o $@ -lm

test: constraint
	./constraint

clean:
	rm -f constraint
