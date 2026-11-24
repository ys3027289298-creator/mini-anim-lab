CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: curve

curve: src/curve.c
	$(CC) $(CFLAGS) $< -o $@ -lm

swing: src/swing.c src/curve.c
	$(CC) $(CFLAGS) src/swing.c -o $@ -lm

test: curve
	./curve

analysis: swing
	./swing

clean:
	rm -f curve swing
