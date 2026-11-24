ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)

.PHONY: all test clean

all: test

# Same entry point as ./run.sh; kept so the whole repo answers to `make test`.
test:
	$(ROOT)/selfcheck/run.sh

clean:
	@:
