using GraphMatFun,BenchmarkTools;
include("setup_graphs.jl");


for (i,g)=enumerate(graphs)
    n=names[i];
    println("# Generating c-files for $n");
    println("# (Ubuntu installation of LAPACKE: sudo apt install liblapacke-dev)");
    fname_mkl=string(tempdir(),"/$(n)_MKL");
    fname_openblas=string(tempdir(),"/$(n)_OpenBLAS");
    gen_code("$(fname_mkl).c",g,funname="$(n)_MKL",lang=GraphMatFun.LangC_MKL(),priohelp=priohelp);
    gen_code("$(fname_openblas).c",g,funname="$(n)_OpenBLAS",lang=GraphMatFun.LangC_OpenBLAS(),priohelp=priohelp);

    println("cat ~/jobb/src/matfun/simulations/cputime/c_main_template.c | sed 's/FUNNAME/d$(n)_MKL/' >> $(fname_mkl).c")

    println("gcc -Wall -lm  -o $(fname_mkl) $(fname_mkl).c -lmkl_rt");
    println("cat ~/jobb/src/matfun/simulations/cputime/c_main_template.c | sed 's/FUNNAME/d$(n)_OpenBLAS/' >> $(fname_openblas).c")
    println("gcc -Wall -lm  -o $(fname_openblas) $(fname_openblas).c -lblas -llapacke");
end
println("");
println("# Run the timing as follows")

for (i,n)=enumerate(names)
    fname_mkl=string(tempdir(),"/$(n)_MKL");
    fname_openblas=string(tempdir(),"/$(n)_OpenBLAS");

    println("$fname_mkl");
    println("echo -n sleeping; sleep 2; echo n"); # Avoid overheating
    println("$fname_openblas");
    println("echo -n sleeping; sleep 2; echo n"); # Avoid overheating
end
