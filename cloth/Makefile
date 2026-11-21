CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pthread

all: cloth

cloth: src/cloth.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: cloth
	./cloth

clean:
	rm -f cloth
