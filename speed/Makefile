CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: speed

speed: src/speed.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: speed
	./speed

clean:
	rm -f speed
