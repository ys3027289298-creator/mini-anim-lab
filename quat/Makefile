CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: quat

quat: src/quat.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: quat
	./quat

clean:
	rm -f quat
