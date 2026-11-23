CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
COMMON := ../../common
ROOT := ../..

MOD_SRCS := \
	$(ROOT)/skeleton/src/skeleton.c \
	$(ROOT)/ik/src/ik.c \
	$(ROOT)/curve/src/curve.c \
	$(ROOT)/constraint/src/constraint.c

BINS := bench_skeleton bench_ik bench_curve bench_constraint

all: $(BINS)

define bench_rule
bench_$(1): bench_$(1).c $(MOD_SRCS) $(COMMON)/src/solver.c
	$$(CC) $$(CFLAGS) -I$(COMMON)/include -DMAL_NO_MAIN -DMAL_UNIFIED=1 \
		bench_$(1).c $(MOD_SRCS) $(COMMON)/src/solver.c -o $$@ -lm
endef
$(eval $(call bench_rule,skeleton))
$(eval $(call bench_rule,ik))
$(eval $(call bench_rule,curve))
$(eval $(call bench_rule,constraint))

test: all
	@for b in $(BINS); do echo "===== $$b ====="; ./$$b || exit 1; done

clean:
	rm -f $(BINS)
