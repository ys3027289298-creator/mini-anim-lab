CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: capstan

capstan: src/capstan.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: capstan
	./capstan

clean:
	rm -f capstan
