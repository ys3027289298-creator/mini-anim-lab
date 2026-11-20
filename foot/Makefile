CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: foot

foot: src/foot.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: foot
	./foot

clean:
	rm -f foot
