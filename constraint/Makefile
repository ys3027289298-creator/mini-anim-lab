CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint

constraint: src/constraint.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: constraint
	./constraint

clean:
	rm -rf constraint out
