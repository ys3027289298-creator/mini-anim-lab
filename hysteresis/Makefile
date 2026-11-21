CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hysteresis

hysteresis: src/hysteresis.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hysteresis
	./hysteresis

clean:
	rm -f hysteresis
