CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: redik

redik: src/redik.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: redik
	./redik

clean:
	rm -f redik
