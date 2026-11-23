CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lighting

lighting: src/lighting.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lighting
	mkdir -p out
	./lighting --out out
	./lighting --no-blend --out out

clean:
	rm -f lighting
	rm -rf out
