CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: camera

camera: src/camera.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: camera
	mkdir -p out
	./camera

clean:
	rm -f camera
	rm -rf out
