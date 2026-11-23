CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
MAL_UNIFIED ?= 1
COMMON := ../common

all: skeleton

skeleton: src/skeleton.c $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include -DMAL_UNIFIED=$(MAL_UNIFIED) \
		src/skeleton.c $(COMMON)/src/solver.c -o $@ -lm

test: skeleton
	./skeleton

clean:
	rm -f skeleton
