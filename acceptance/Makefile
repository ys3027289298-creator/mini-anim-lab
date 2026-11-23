CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
COMMON := ../common
ROOT := ..

MOD_SRCS := \
	$(ROOT)/skeleton/src/skeleton.c \
	$(ROOT)/ik/src/ik.c \
	$(ROOT)/curve/src/curve.c \
	$(ROOT)/constraint/src/constraint.c

all: invalid

invalid: src/invalid.c $(MOD_SRCS) $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include \
		-DMAL_NO_MAIN -DMAL_UNIFIED=1 \
		src/invalid.c $(MOD_SRCS) $(COMMON)/src/solver.c -o $@ -lm

test: invalid
	./invalid

clean:
	rm -f invalid
