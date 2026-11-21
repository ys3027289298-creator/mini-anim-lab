CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

CHAIN_SRC = ../constraint/src/chain.c
CHAIN_HDR = ../constraint/src/chain.h

all: audit

audit: src/audit.c $(CHAIN_SRC) $(CHAIN_HDR)
	$(CC) $(CFLAGS) src/audit.c $(CHAIN_SRC) -I../constraint/src -o $@ -lm

test: audit
	./audit

clean:
	rm -f audit
