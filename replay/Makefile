CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: replay

replay: src/replay.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: replay
	@tmp=$$(mktemp -d); \
	sha256sum samples/* tests/* > $$tmp/before.sha; \
	echo "--- check: determinism + reference bound ---"; \
	./replay check samples/chain3.anim samples/chain3.ref; \
	echo "--- cross-run: record twice, diff ---"; \
	./replay record samples/chain3.anim $$tmp/a.rec; \
	./replay record samples/chain3.anim $$tmp/b.rec; \
	./replay diff $$tmp/a.rec $$tmp/b.rec; \
	echo "--- invalid inputs must be rejected ---"; \
	! ./replay check tests/bad_nan.anim samples/chain3.ref 2>$$tmp/err1; grep -q "NaN" $$tmp/err1 && echo "nan_rejected=1"; \
	! ./replay check tests/bad_negt.anim samples/chain3.ref 2>$$tmp/err2; grep -q "negative" $$tmp/err2 && echo "negt_rejected=1"; \
	echo "--- samples untouched ---"; \
	sha256sum -c $$tmp/before.sha; \
	rm -rf $$tmp

clean:
	rm -f replay
