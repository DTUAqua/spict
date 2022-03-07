R?=R

all:
	@if [ -f "ignore" ]; then echo -n $(shell basename $(shell pwd)); echo "... ignored"; else $(MAKE) doit; fi
doit:
	rm -f res.out
	unset MAKEFLAGS; echo 'source("script.R")' | ${R} --slave > /dev/null
	@echo -n $(shell basename $(shell pwd))
	@echo -n "... "
	@diff -q --ignore-all-space res.out res.EXP; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
	  echo "OK"; \
	else \
	  echo "Error:"; \
	  diff --ignore-all-space res.out res.EXP;\
	fi
