using GraphMatFun
using GenericSVD,LinearAlgebra

function get_taylor(eltype,n)
    c=1 ./factorial.(convert.(big(real(eltype)),0:(n-1)));
    c=convert.(eltype,c);
end
function get_lsqr(n,discr)
    Z=zeros(big(eltype(discr)),size(discr,1),n);
    Z[:,1]=ones(size(discr,1));
    for s=2:n
        Z[:,s]=(Z[:,s-1]).*(discr);
    end
    c=Z\exp.(big.(discr));
    return convert.(eltype(discr),c);
end


mutable struct Simulation
    m # Matrix-Matrix multiplications
    n # Nof samples
    f # function
    rho # disc radius
    init  # way to initialize # :lsqr or :taylor or :prev
    graph # :mono, ..
    eltype # ComplexF64
    opt_kwargs # Dict with override kwargs for optimization
end

function Simulation(m;n=100,f=exp,rho=1.0,init=:prev,graph=:prev,eltype=ComplexF64,opt_kwargs=Dict())
    return Simulation(m,n,f,rho,init,graph,eltype,opt_kwargs);
end
function initsim(s,prev_graph,prev_cref)
    n=s.n
    discr=s.rho*exp.(1im*range(0,2*pi,length=s.n)[1:end-1]);

    discr=convert.(s.eltype,discr);

    if (s.f != exp)
        error("Function not implemented");
    end


    c=Vector();
    graph=Nothing;
    cref=Nothing;
    if (s.graph == :prev)
        graph=Compgraph(s.eltype, prev_graph);
        cref=prev_cref;
    else
        if (s.graph == :sid) # generator has nof mult as input
                (graph,cref)=gen_sid_exp(m)
        else
            deg0=m-1;
            graph_mult=0;
            oldgraph=Nothing
            oldcref=Nothing;
            while graph_mult<m+1

                @show deg0
                deg0 += 1;
                if (s.init==:lsqr)
                    bdiscr=big.(discr);
                    c=get_lsqr(deg0+1,bdiscr)
                    c=convert.(s.eltype,c);
                elseif (s.init==:taylor)
                    c=get_taylor(s.eltype,deg0+1)
                elseif (s.init == :prev)
                else
                    error("Unknown init $(s.init)")

                end
                oldgraph=graph;
                oldcref=cref;
                if (s.graph == :mono)
                    (graph,cref)=gen_monomial_degopt(c)
                elseif (s.graph == :horner)
                    (graph,cref)=gen_horner_degopt(c)
                elseif (s.graph == :ps)
                    (graph,cref)=gen_ps_degopt(c)
#                elseif (s.graph == :sid)
#                    (graph,cref)=gen_sid_exp(m)
                else
                    error("Wrong graph type $(s.graph)");
                end
                graph_mult=count(values(graph.operations) .==:mult);

                @show graph_mult
            end

            graph=oldgraph; # Useful at termination
            cref=oldcref;


            println("Using degree $(deg0-1) graph: $(typeof(graph))");
        end

    end



    if (m != count(values(graph.operations) .==:mult))
        @show m, count(values(graph.operations) .==:mult)
        error("Incorrect number of multiplications")
    end

    graph=Compgraph(s.eltype,graph);
    return (graph,cref,discr);
end

function runsim(s::Simulation,prev_graph,prev_cref)
    (graph,cref,discr)=initsim(s,prev_graph,prev_cref)
    println("Running graph \"$(s.graph)\" with $(s.m) multiplications rho=$(s.rho) ");

    opt_gauss_newton!(graph,exp,discr;logger=1,
                      stoptol=1e-16,cref=cref,
                      s.opt_kwargs...);
    return (deepcopy(graph),cref);
end


function showerr(s,graph)
    n=s.n
    rho=s.rho
    discr=s.rho*exp.(1im*range(0,2*pi,length=s.n)[1:end-1]);

    discr=convert.(s.eltype,discr);

    if (s.f != exp)
        error("Function not implemented");
    end

    err=norm((s.f.(discr)-eval_graph(graph,discr))./s.f.(discr),Inf)
    println("Target error: $err");
    return err
end



function make_accurate(sim,sim_target)
    new_sim=deepcopy(sim);
    new_sim.n=sim_target.n;
    return new_sim;
end

function make_agressive(sim)
    new_sim=deepcopy(sim);
    new_sim.opt_kwargs[:droptol]=new_sim.opt_kwargs[:droptol]/10;

    return new_sim;
end
