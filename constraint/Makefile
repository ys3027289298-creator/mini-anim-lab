CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: constraint topo

constraint: src/constraint.c
	$(CC) $(CFLAGS) $< -o $@ -lm

topo: src/topo.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: constraint topo
	./constraint
	./topo

clean:
	rm -f constraint topo
