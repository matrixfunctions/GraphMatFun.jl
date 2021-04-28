SHELL := /bin/bash
.SILENT: all create_dir generate \
	runall run_c run_julia run_matlab \
	cleanall clean_c clean_julia clean_matlab clean_output

N=1000 # Default matrix size
TEMPDIR=$(realpath ../generated_code/)

# Clean the /tmp/ directory, then generate the code, and run the experiments.
all: cleanall generate runall




# Create temporary directory if not already existing.
create_dir:
	if [ ! -d $(TEMPDIR) ]; then \
		mkdir $(TEMPDIR); \
	fi

# Generate and compile test files in the temporary directory.
# Nothing happens in the current directory.
generate: create_dir
	echo "Generating C, Julia, and Matlab code..."
	export TMPDIR="$(TEMPDIR)"; cd ../..; \
		julia simulations/cputime/gen_cputime_all_lang.jl  &> \
		$(TEMPDIR)/gen_all_lang.jl.output
	echo "Compiling generated C code..."
	# For MKL as part of Intel OneAPI, setting the environment correctly
	# requires running
	#     source /opt/intel/oneapi/setvars.sh
	# from the parent shell.
	gcc -Wall -lm -o $(TEMPDIR)/run_cputime_c_MKL \
		$(TEMPDIR)/run_cputime_c_MKL.c -lmkl_rt
	gcc -Wall -lm  -o $(TEMPDIR)/run_cputime_c_OpenBLAS \
		$(TEMPDIR)/run_cputime_c_OpenBLAS.c -lblas -llapacke




# Run experiments generated in the temporary directory.
runall: run_c run_julia run_matlab

run_c:
	if [ `cat /sys/devices/system/cpu/intel_pstate/no_turbo`!="1" ]; then \
		echo "Warning: CPU-time boosting seems to be enabled."; \
		echo "Disable on ubuntu: https://askubuntu.com/a/620114"; \
	fi
	$(TEMPDIR)/run_cputime_c_MKL &> \
		$(TEMPDIR)/run_cputime_c_MKL.output
	$(TEMPDIR)/run_cputime_c_OpenBLAS &> \
		$(TEMPDIR)/run_cputime_c_OpenBLAS.output

run_julia:
	cd $(TEMPDIR)/; julia run_cputime_julia_MKL.jl &> \
		$(TEMPDIR)/run_cputime_julia_MKL.jl.output
	cd $(TEMPDIR)/; julia run_cputime_julia_OpenBLAS.jl &> \
		$(TEMPDIR)/run_cputime_julia_OpenBLAS.jl.output

run_matlab:
	cd $(TEMPDIR)/; matlab -batch "run_cputime_matlab_OpenBLAS" &> \
		$(TEMPDIR)/run_cputime_matlab_OpenBLAS.m.output




# Clean the temporary directory.
cleanall: clean_c clean_julia clean_matlab clean_output

clean_c:
	echo "Removing C source and executable files..."
	rm -f $(TEMPDIR)/*.c $(TEMPDIR)/run_cputime_c_*

clean_julia:
	echo "Removing Julia files..."
	rm -f $(TEMPDIR)/*.jl

clean_matlab:
	echo "Removing Matlab files..."
	rm -f $(TEMPDIR)/*.m

clean_output:
	echo "Removing output files"
	rm -f $(TEMPDIR)/*.output
