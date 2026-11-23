CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
MAL_UNIFIED ?= 1
COMMON := ../common

all: curve

curve: src/curve.c $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include -DMAL_UNIFIED=$(MAL_UNIFIED) \
		src/curve.c $(COMMON)/src/solver.c -o $@ -lm

test: curve
	./curve

clean:
	rm -f curve
