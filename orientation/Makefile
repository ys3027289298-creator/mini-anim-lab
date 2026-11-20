CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: orientation

orientation: src/orientation.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: orientation
	./orientation

clean:
	rm -f orientation
