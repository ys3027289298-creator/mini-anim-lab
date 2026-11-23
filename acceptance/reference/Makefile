CC ?= gcc
# -Wno-error only; the original ik.c predates its <string.h> include and
# emits a warning under newer gcc. It builds and runs correctly, and the
# point of this target is to reproduce the original bit-for-bit behavior.
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wno-implicit-function-declaration \
          -Wno-builtin-declaration-mismatch

REFS := ref_skeleton ref_ik ref_curve ref_constraint

all: $(REFS)

ref_%: ref_%.c
	$(CC) $(CFLAGS) -I. $< -o $@ -lm

test: all
	@for r in $(REFS); do echo "===== $$r ====="; ./$$r || exit 1; done

clean:
	rm -f $(REFS)
