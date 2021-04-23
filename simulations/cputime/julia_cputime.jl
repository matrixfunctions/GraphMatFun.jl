using GraphMatFun,BenchmarkTools,LinearAlgebra;
include("setup_graphs.jl");

@show names
function get_cost(graph::Compgraph)
    return count(values(graph.operations) .== :mult)+
           (4/3)*count(values(graph.operations) .== :ldiv)

end

println("Matrix norm: $(opnorm(A0,1))");
if (!isdefined(BLAS,:get_config))
    println("BLAS vendor: ",BLAS.vendor())
else
    println("BLAS config: ",BLAS.get_config())
end
graphs=Vector{Any}(graphs);
pushfirst!(graphs,:exp);
pushfirst!(names,"exp");

for (i,g)=enumerate(graphs)
    local n;

    n=names[i];
    A=deepcopy(A0);

    print("$n ");
    if (g isa Compgraph)
        print("cost $(get_cost(g)) ");

        gen_code(string(tempdir(),"/exp_$n.jl"),g,funname="exp_$n",priohelp=priohelp);
        gen_code(string(tempdir(),"/exp_$n.m"),g,funname="exp_$n",lang=LangMatlab(),priohelp=priohelp);

        include(string(tempdir(),"/exp_$n.jl"));

        print("time: ");
        s=Symbol("exp_$n");
        bb=@benchmark eval(Expr(:call,$s,$A))
    else
        print("time: ");
        bb=@benchmark exp($A);
    end


    #@show bb.times*1e-9
    mm=median(bb.times)*1e-9;
    println("$(mm) mem: $(bb.memory)");
    print("      ")
    @show bb.times*1e-9

    #println("median: mem: $(bb.memory)");
end
