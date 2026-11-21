CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sleep

sleep: src/sleep.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sleep
	./sleep

clean:
	rm -f sleep
