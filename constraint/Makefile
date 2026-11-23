CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
MAL_UNIFIED ?= 1
COMMON := ../common

all: constraint

constraint: src/constraint.c $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include -DMAL_UNIFIED=$(MAL_UNIFIED) \
		src/constraint.c $(COMMON)/src/solver.c -o $@ -lm

test: constraint
	./constraint

clean:
	rm -f constraint
