CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
MAL_UNIFIED ?= 1
COMMON := ../common

all: ik

ik: src/ik.c $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include -DMAL_UNIFIED=$(MAL_UNIFIED) \
		src/ik.c $(COMMON)/src/solver.c -o $@ -lm

test: ik
	./ik

clean:
	rm -f ik
