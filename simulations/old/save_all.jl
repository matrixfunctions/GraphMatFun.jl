include("set_lists_all.jl");

if (isnothing(target))
    error("Target not set")
end

println("Saving simulation");

rho=target.rho;
println("Domain: $(rho), m=",target.m);
for (i,g)=enumerate(graph_list)
    name=name_list[i];
    if (!isnothing(g))
        rho_str=replace(string(rho),"."=>"_");
        fname=joinpath("simulations","graphs","exp_m$(m)_$(name)_$(rho_str).cgr");
        err=Float64(showerr(target,g,false))
        get_coeffs(g)

        isreal=false;
        if (norm(imag.(get_coeffs(g)))==0)
            isreal=true;
            g=Compgraph(real(eltype(g)),g);
        end
        println("Saving $name $fname isreal=$isreal");
        export_compgraph(g, fname;
                         fun="exp", dom="$(rho)D", err="$err",
                         order=get_topo_order(g)[1],
                         descr="Matrix exponential with degree optimal form optimized starting with $name");

    else
        println("Skip saving $name");
    end
end
