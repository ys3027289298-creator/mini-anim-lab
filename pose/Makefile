CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: pose

pose: src/pose.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: pose
	./pose

clean:
	rm -f pose
