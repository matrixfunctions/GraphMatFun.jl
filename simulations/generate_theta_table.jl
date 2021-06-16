using LinearAlgebra,GraphMatFun,Printf


mutable struct ThetaEntry
    mul
    rstr
    theta_val # starting value for fzero and result
    method
    params
    err
    rho
end

function compute_err(graph,rho,n=400)

    f=x->exp(x)
    nn=n;
    discr=rho*exp.(1im*range(0,2*pi,length=nn)[1:end-1]);
    discr=convert.(Complex{BigFloat},discr);
    err=norm((f.(discr)-eval_graph(graph,discr))./f.(discr),Inf)
    return err;
end
function latexname(str)
    replace(str,"+"=>"&");
end


mv=[4; 5; 5; 6; 6; 7; 7; 8]
rv=["0_69"; "1_68"; "1_9"; "2_22"; "2_7"; "3_59"; "6_0"; "13_5"];
methods=["mono+GN"; "ps+GN"; "sastre+GN"; "bbc+GN"; "sid+GN"; "sastre"; "bbc"; "sid" ];

simulations=Matrix{ThetaEntry}(undef,size(methods,1),size(mv,1));
for (i,m)=enumerate(methods)
    for (j,r)=enumerate(rv)
        mul=mv[j];
        theta_val=parse(Float64,replace(r,"_" => "."));
        simulations[i,j]=ThetaEntry(mul,r,theta_val/1.5,m,Dict(),0,theta_val);

    end
end
#
simulations[1,8].theta_val=14;
simulations[end,8].theta_val=5;
simulations[3,6].theta_val=7.0;
#simulations[3,3].theta_val=0.5;
skip_theta=false;
println("Theta table");
for i=1:size(simulations,1);
    method=simulations[i,1].method
    print(" $(latexname(method)) ");
    for j=1:size(simulations,2);

        sim=simulations[i,j];
        mul=sim.mul


        theta_val=sim.theta_val;
        rstr=sim.rstr;
#        if (mul<7) # Hack to only print high
#            continue;
#        end;

        basename="simulations/newgraphs/exp";
        fname="$(basename)_m$(mul)_$(method)_$(rstr).cgr";

        if (!isfile(fname))
            print(" & X");
            continue
        end
        graph=import_compgraph(fname);
        sim.err=compute_err(graph,sim.rho);
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
            if (!skip_theta)
            theta=compute_bwd_theta_exponential(graph,coefftype=BigFloat,
                                                tolerance=eps()/2,theta_init=big(theta_val))[2];
            end


        catch e
        end

        sim.theta_val=theta;
        if (!isnan(theta))
            thetastr=@sprintf("%.3f",theta);
        else
            thetastr="x";
        end


        print(" & $thetastr")

    end
    println("  \\\\");
end
println("Error table");
for i=1:size(simulations,1);
    method=simulations[i,1].method
    print(" $(latexname(method)) ");
    for j=1:size(simulations,2);
        sim=simulations[i,j];
        e=sim.err
        if (e>1 || e == 0)
            estr = "\$\\times\$"
        else
            estr=@sprintf("%.1e",e);
        end

        print(" & $estr");


    end
    println("\\\\");

end



#
#for (i,m)=enumerate(methods)
#    print(" $m ");
#    for (j,r)=enumerate(rv)
#        mul=mv[j];
#        if (mul<8) # Hack to only print high
#            continue;
#        end;
#
#
#        rval=parse(Float64,replace(r,"_" => "."));
#
#
#        rval=12.5;
#
#        basename="simulations/graphs/exp";
#        fname="$(basename)_m$(mul)_$(m)_$(r).cgr";
#        if (!isfile(fname))
#            print("X");
#            continue
#        end
#        graph=import_compgraph(fname);
#        x=get_coeffs(graph);
#        if (norm(imag.(x))>eps())
#            warning("Non-real")
#        else
#            x=real.(x);
#            set_coeffs!(graph,x);
#            graph=Compgraph(BigFloat,graph)
#        end
#        theta=NaN
#        try
#            # If the fzero function does not find a root it
#            # will throw an error
#            theta=compute_bwd_theta_exponential(graph,coefftype=BigFloat,tolerance=eps()/2,theta_init=big(rval))[2];
#        catch e
#            throw(e)
#        end
#
#        thetastr=@sprintf("%.3f",theta);
#        print(" & $thetastr")
#    end
#    println("  \\\\");
#end
#
