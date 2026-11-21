CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint equiv

constraint: src/constraint.c
	$(CC) $(CFLAGS) $< -o $@ -lm

equiv: src/equiv.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: constraint equiv
	./constraint
	./equiv
	! ./equiv 0

clean:
	rm -f constraint equiv
