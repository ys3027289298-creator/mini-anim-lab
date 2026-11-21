CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: wave

wave: src/wave.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: wave
	./wave

clean:
	rm -f wave
