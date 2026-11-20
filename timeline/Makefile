CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: timeline

timeline: src/timeline.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: timeline
	./timeline

clean:
	rm -f timeline
