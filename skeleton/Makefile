CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: skeleton

skeleton: src/skeleton.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: skeleton
	./skeleton

clean:
	rm -f skeleton
