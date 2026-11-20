CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: grasp

grasp: src/grasp.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: grasp
	./grasp

clean:
	rm -f grasp
