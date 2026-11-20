CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: palette

palette: src/palette.c
	$(CC) $(CFLAGS) $< -o $@ -lm

palette-thread: src/palette.c
	$(CC) $(CFLAGS) -DPAL_WITH_PTHREAD -pthread $< -o $@ -lm

# On WSL2 / newer kernels TSan needs the legacy address-space layout.
# setarch is harmless where not needed; suppressions cover the seqlock
# payload overlap that TSan cannot model (see tsan-suppressions.txt).
palette-tsan: src/palette.c tsan-suppressions.txt
	$(CC) -O1 -g $(CFLAGS) -fsanitize=thread -DPAL_WITH_PTHREAD -pthread $< -o $@ -lm

test: palette
	./palette

thread-test: palette-thread
	./palette-thread --thread

tsan-test: palette-tsan
	TSAN_OPTIONS="halt_on_error=0 suppressions=tsan-suppressions.txt" \
	  setarch $(shell uname -m) --addr-no-randomize ./palette-tsan --thread

clean:
	rm -f palette palette-thread palette-tsan

.PHONY: all test thread-test tsan-test clean
