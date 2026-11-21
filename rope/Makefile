CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: rope

rope: src/rope.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: rope
	./rope

clean:
	rm -f rope
