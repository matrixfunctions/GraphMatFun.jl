using GraphMatFun;
include("setup_graphs.jl");




function gen_julia_main(N,col,graphs,names,mkl)
    fname="simulations/cputime/julia_main_template.jl";
    lines=readlines(fname);
    start_repeated_line=0;
    end_repeated_line=0;
    for (i,line)=enumerate(lines)
        line=replace(line, "COLUMN" => string(col));
        if (mkl)
            line=replace(line, "USING" => "using MKL;");
        else
            line=replace(line, "USING" => "# NO MKL");
        end

        lines[i]=line
        if contains(line,"START REPEATED CODE")
            start_repeated_line=i;
        end
        if contains(line,"END REPEATED CODE")
            end_repeated_line=i;
        end

    end


    repeated_code=lines[start_repeated_line+1:end_repeated_line-1];

    postprocess=lines[(end_repeated_line+1):end];
    lines=lines[1:(start_repeated_line-1)]
    for (i,g)=enumerate(graphs)
        local n;
        if (g isa Compgraph || (g == :exp_julia))
            n=names[i];
            this_sim_code=deepcopy(repeated_code);
            for (j,line)=enumerate(this_sim_code)
                if (g isa Compgraph)
                    line=replace(line,"INCLUDE" =>
                      string("include(joinpath(tempdir(),\"NAME.jl\"));"));

                else
                    line=replace(line,"INCLUDE" => "exp_julia=exp");
                end

                line=replace(line,"NAME" => n);
                line=replace(line,"FUNCTION" => n);

                this_sim_code[j]=line;
            end

            # Append
            map(code->push!(lines,code), this_sim_code)
        else
            push!(lines,"# SKIPPING $g");
            push!(lines,"");
        end


    end


    lines=[lines;postprocess]


    return lines;

end
function get_cost(graph::Compgraph)
    return count(values(graph.operations) .== :mult)+
           (4/3)*count(values(graph.operations) .== :ldiv)

end



println("Generating graphs");
for (i,g)=enumerate(graphs)
    local n;

    n=names[i];

    print("Generating code for $n ");
    if (g isa Compgraph)
        print("cost $(get_cost(g)) ");
        gen_code(string(tempdir(),"/$n.jl"),g,funname="$n",priohelp=priohelp);
        gen_code(string(tempdir(),"/$n.m"),g,funname="$n",lang=LangMatlab(),priohelp=priohelp);

        # Distinguish between MKL and OpenBLAS in the graph generators
        # already for C. Since includes etc are different.
        fname_mkl=string(tempdir(),"/$(n)_MKL");
        fname_openblas=string(tempdir(),"/$(n)_OpenBLAS");
        gen_code("$(fname_mkl).c",g,funname="$(n)_MKL",lang=GraphMatFun.LangC_MKL(),priohelp=priohelp);
        gen_code("$(fname_openblas).c",g,funname="$(n)_OpenBLAS",lang=GraphMatFun.LangC_OpenBLAS(),priohelp=priohelp);

    end
    println();
end

lines=gen_julia_main(2000,1,graphs,names,true);
open(string(tempdir(),"/run_cputime_julia_MKL.jl"),"w") do io
    for line=lines
        write(io,line,"\n");
    end
end

lines=gen_julia_main(2000,1,graphs,names,false);
open(string(tempdir(),"/run_cputime_julia_OpenBLAS.jl"),"w") do io
    for line=lines
        write(io,line,"\n");
    end
end
