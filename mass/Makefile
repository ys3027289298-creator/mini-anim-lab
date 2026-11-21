CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: mass

mass: src/mass.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: mass
	./mass

clean:
	rm -f mass
