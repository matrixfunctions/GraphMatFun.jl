using MAT,GraphMatFun,BenchmarkTools,LinearAlgebra;
include("setup_graphs.jl");

graph_m6=graphs[1];
graph_m7=graphs[2];


gen_code("/tmp/exp_m6.jl",graph_m6,funname="exp_m6",priohelp=priohelp);
gen_code("/tmp/exp_m6.m",graph_m6,funname="exp_m6",lang=LangMatlab());

gen_code("/tmp/exp_m7.jl",graph_m7,funname="exp_m7");
gen_code("/tmp/exp_m7.m",graph_m7,funname="exp_m7",lang=LangMatlab());



include("/tmp/exp_m6.jl");
include("/tmp/exp_m7.jl");
println("Timing");

println("timing m7");
A=deepcopy(A0)
@btime exp_m7(A)
println("timing m6");
A=deepcopy(A0)
@btime exp_m6(A)
A=deepcopy(A0);
@btime exp(A);
