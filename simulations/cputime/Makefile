SHELL := /bin/bash

all:
	cd ../..; julia simulations/cputime/c_cputime.jl
	cp c_main_template.c /tmp/
	for file in `ls -d /tmp/* | grep -P 'MKL|OpenBLAS'`; do \
		filebasename=`basename $$file`; \
		filename=$${filebasename%.*}; \
		cat /tmp/c_main_template.c | \
			sed "s/FUNNAME/d$$filename/" \
			>> $$file; \
		type="$${filename#"$${filename%_*}_"}"; \
		if [ $$type == "MKL" ]; then \
			gcc -Wall -lm -o /tmp/$$filename $$file -lmkl_rt; \
		else \
			gcc -Wall -lm  -o /tmp/$$filename $$file -lblas -llapacke; \
		fi; \
	done
