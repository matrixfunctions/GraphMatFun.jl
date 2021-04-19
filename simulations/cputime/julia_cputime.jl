using MAT,GraphMatFun,BenchmarkTools,LinearAlgebra;
include("setup_graphs.jl");

graph_m6=graphs[1];
graph_m7=graphs[2];


gen_code("/tmp/exp_m6.jl",graph_m6,funname="exp_m6",priohelp=priohelp);
gen_code("/tmp/exp_m6.m",graph_m6,funname="exp_m6",lang=LangMatlab(),priohelp=priohelp);

gen_code("/tmp/exp_m7.jl",graph_m7,funname="exp_m7",priohelp=priohelp);
gen_code("/tmp/exp_m7.m",graph_m7,funname="exp_m7",lang=LangMatlab(),priohelp=priohelp);

println("BLAS: ",BLAS.vendor())
for (i,g)=enumerate(graphs)
    local n;
    n=names[i];

    println("Generating $n");

    gen_code("/tmp/exp_$n.jl",g,funname="exp_$n",priohelp=priohelp);
    gen_code("/tmp/exp_$n.m",g,funname="exp_$n",lang=LangMatlab(),priohelp=priohelp);

    include("/tmp/exp_$n.jl");

    println("Timing $n");
    s=Symbol("exp_$n");
    A=deepcopy(A0);
    bb=@benchmark eval(Expr(:call,$s,:A0))
    println("median: $(bb.times[2]*1e-9) mem: $(bb.memory)");
end
