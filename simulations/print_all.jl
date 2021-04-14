
include("set_lists_all.jl");

if (isnothing(target))
    error("Target not set")
end

println("Printing simulation");

println("Domain: $(target.rho), m=",target.m);
for (i,g)=enumerate(graph_list)
    name=name_list[i];
    if (!isnothing(g))
        err=showerr(target,g,false)
        println("$name $(Float64(err))");
    else
        println("$name x");
    end
end
