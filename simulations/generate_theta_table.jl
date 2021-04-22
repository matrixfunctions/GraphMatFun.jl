using LinearAlgebra,GraphMatFun,Printf


mv=[4; 5; 5; 6; 6; 7; 7]
rv=["0_69"; "1_68"; "1_9"; "2_22"; "2_7"; "3_59"; "6_0"];
methods=["mono_taylor"; "PS_taylor"; "Sastre+"; "BBC+"; "SID+"; "Sastre"; "BBC"; "SID" ];

for (i,m)=enumerate(methods)
    print(" $m ");
    for (j,r)=enumerate(rv)
        mul=mv[j];
        rval=parse(Float64,replace(r,"_" => "."));

        basename="simulations/graphs/exp";
        fname="$(basename)_m$(mul)_$(m)_$(r).cgr";
        if (!isfile(fname))
            print("X");
            continue
        end
        graph=import_compgraph(fname);
        x=get_coeffs(graph);
        if (norm(imag.(x))>eps())
            warning("Non-real")
        else
            x=real.(x);
            set_coeffs!(graph,x);
            graph=Compgraph(BigFloat,graph)
        end
        theta=NaN
        try
            # If the fzero function does not find a root it
            # will throw an error
            theta=compute_bwd_theta_exponential(graph,coefftype=BigFloat,tolerance=eps()/2,theta_init=big(rval))[2];
        catch
        end

        thetastr=@sprintf("%.3f",theta);
        print(" & $thetastr")
    end
    println("  \\\\");
end
