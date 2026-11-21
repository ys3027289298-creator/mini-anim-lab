CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: frame

frame: src/frame.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: frame
	./frame

clean:
	rm -f frame
