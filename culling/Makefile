CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: culling

culling: src/culling.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: culling
	./culling

clean:
	rm -f culling
