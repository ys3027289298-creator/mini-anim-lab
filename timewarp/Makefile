CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: timewarp

timewarp: src/timewarp.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: timewarp
	./timewarp

clean:
	rm -f timewarp
