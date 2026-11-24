CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lod

lod: src/lod.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# Run twice and diff: the whole report, trigger frames included, must be
# byte-identical between processes.
test: lod
	@./lod > .lod_run1; s1=$$?; ./lod > .lod_run2; s2=$$?; \
	if cmp -s .lod_run1 .lod_run2; then echo repro_binary=1; \
	else echo repro_binary=0; fi; \
	cat .lod_run1; rm -f .lod_run1 .lod_run2; \
	test $$s1 -eq 0 -a $$s2 -eq 0

clean:
	rm -f lod .lod_run1 .lod_run2
