CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: curve

curve: src/curve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: curve
	./curve

clean:
	rm -f curve
