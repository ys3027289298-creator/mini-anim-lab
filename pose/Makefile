CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
COMMON := ../common
ROOT := ..

MOD_SRCS := \
	$(ROOT)/skeleton/src/skeleton.c \
	$(ROOT)/ik/src/ik.c \
	$(ROOT)/curve/src/curve.c \
	$(ROOT)/constraint/src/constraint.c

all: pose

pose: src/pose.c $(MOD_SRCS) $(COMMON)/src/solver.c
	$(CC) $(CFLAGS) -I$(COMMON)/include \
		-DMAL_NO_MAIN -DMAL_UNIFIED=1 \
		src/pose.c $(MOD_SRCS) $(COMMON)/src/solver.c -o $@ -lm

test: pose
	./pose

clean:
	rm -f pose
