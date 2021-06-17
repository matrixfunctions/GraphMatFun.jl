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


    if (s.graph isa Tuple)
        # Starting graph is hardcoded
        return (s.graph[1],s.graph[2],discr)
    end

    c=Vector();
    graph=Nothing;
    cref=Nothing;
    if (s.graph == :prev)
        graph=Compgraph(s.eltype, prev_graph);
        cref=prev_cref;
    else
        if (s.graph in [:sid,:bbc,:sastre]) # generator has nof mult as input

            if (s.graph == :sid)
                (graph,cref)=gen_sid_exp(m)
            elseif (s.graph == :bbc)
                (graph,cref)=gen_bbc_basic_exp(m)
            elseif (s.graph == :sastre)
                (graph,cref)=GraphMatFun.gen_sastre_basic_exp(m,:y1s)
            end
        else
            deg0=m-1;
            @show deg0
            graph_mult=0;
            oldgraph=Nothing
            oldcref=Nothing;
            # Automatically determine the order
            while graph_mult<m+1

                deg0 += 1;
                if (s.init==:lsqr)
                    bdiscr=big.(discr);
                    c=get_lsqr(deg0+1,bdiscr)
                    c=convert.(s.eltype,c);
                elseif (s.init==:taylor)
                    c=get_taylor(s.eltype,deg0)
                elseif (s.init == :prev)
                else
                    error("Unknown init $(s.init)")

                end
                oldgraph=graph;
                oldcref=cref;
                if (s.graph == :mono)
                    (graph,cref)=gen_monomial_degopt(c)
                elseif (s.graph == :mono2)
                    (graph_org,cref_org)=gen_monomial_degopt(c[1:end-1])
                    (graph,cref)=gen_degopt_by_squaring(graph_org);
                elseif (s.graph == :horner)
                    (graph,cref)=gen_horner_degopt(c)
                elseif (s.graph == :ps)
                    (graph,cref)=gen_ps_degopt(c)
                else
                    error("Wrong graph type $(s.graph)");
                end
                graph_mult=count(values(graph.operations) .==:mult);

            end

            # We go one step too far in the automatic determination
            graph=oldgraph;
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
    if (s.graph isa Tuple)
        println("Running graph initialized with pre-computed graph with $(s.m) multiplications rho=$(s.rho) ");

    elseif (s.graph != :prev)
        println("Running graph \"$(s.graph)\" initialized as \"$(s.init)\" with $(s.m) multiplications rho=$(s.rho) ");
    else
        d=s.opt_kwargs[:droptol];
        println("Running with droptol \"$(d)\" with $(s.m) multiplications rho=$(s.rho) ");
    end

    opt_gauss_newton!(graph,exp,discr;logger=1,
                      stoptol=1e-30,cref=cref,
                      s.opt_kwargs...);
    return (deepcopy(graph),cref);
end


function showerr(s,graph,output=true)
    n=s.n
    rho=s.rho
    discr=s.rho*exp.(1im*range(0,2*pi,length=s.n)[1:end-1]);

    discr=convert.(s.eltype,discr);

    if (s.f != exp)
        error("Function not implemented");
    end

    err=norm((s.f.(discr)-eval_graph(graph,discr))./s.f.(discr),Inf)
    imagnorm=Float64(norm(imag.(get_coeffs(graph))));
    if (output)
        if (imagnorm>0)
            imagstr="norm(imag(coeffs))= $imagnorm"
        else
            imagstr="";
        end
        println("Target error: $(Float64(err)) $imagstr");
    end

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


function run_simlist(simlist)
    graph=0; cref=[];
    for sim=simlist
        if (sim isa String)
            println(sim)
        else
            oldgraph=graph
            (graph,cref)=runsim(sim,graph,cref);
            err=showerr(target,graph);
            if (err>1e20)
                println("Too large error. Aborting.");
                break;
            end

        end
    end
    return (graph,cref)
end



function change_param(sim,param,factor)
    sim=deepcopy(sim);
    if (param in [:droptol,:γ0])
        o=sim.opt_kwargs
        o[param]=o[param]*factor;
        println("$param=$(o[param])");
    else
        T=typeof(getfield(sim,param));
        setfield!(sim,param,convert(T,getfield(sim,param)*factor));
        println("$param= $(getfield(sim,param))");
    end

    return sim;
end


function read_one_key(; io = stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io.handle, raw)
    setraw!(true)
    s=read(io, 1)
    setraw!(false)
    return s
end


function interactive_simulations(init_sim,sim0,predefsims="")

    max_j=1000;
    simlist=Vector{Any}(undef,max_j);
    commandlist=Vector{String}(undef,max_j);
    graphlist=Vector{Any}(undef,max_j);
    (graph,cref)=runsim(init_sim,0,0);
    simlist[1]=init_sim;
    graphlist[1]=graph;
    commandlist[1]="I";
    err=showerr(target,graph);
    println("Initialerr:$err");
    (graph,cref)=runsim(sim0,graph,cref);
    simlist[2]=sim0;
    graphlist[2]=graph;
    commandlist[2]="s";
    j=2;
    command="";



    while (command != "q")
        graph=graphlist[j];
        err=showerr(target,graph);
        print("Commandlist : "*String(join(commandlist[1:j]))* " ");

        if (j+1>length(predefsims))
            x=String(read_one_key());
        else
            x=string(predefsims[j+1]);
        end

        println("$x");
        if (x=="b")
            j=j-1;
        else
            skip=false;
            if (x=="s")
                sim=deepcopy(simlist[j]);
            elseif (x=="D")
                sim=change_param(simlist[j],:droptol,2);
            elseif (x=="d")
                sim=change_param(simlist[j],:droptol,1/2);
            elseif (x=="g")
                sim=change_param(simlist[j],:γ0,0.8);
            elseif (x=="G")
                sim=change_param(simlist[j],:γ0,1/0.8);
            elseif (x=="n")
                sim=change_param(simlist[j],:n,1/2);
            elseif (x=="N")
                sim=change_param(simlist[j],:n,2);
            elseif (x=="r")
                graph=deepcopy(graph)
                v=get_coeffs(graph);
                set_coeffs!(graph,real.(v));
                println("Realify (removing $(Float64(norm(imag.(v)))) )");
                sim=deepcopy(simlist[j]);
            elseif (x=="l")
                # Linear squares on degopt coeffs
                graph=deepcopy(graph)
                y_cref=get_degopt_crefs(simlist[j].m)[2];

                yy0=get_coeffs(graph,y_cref);
                opt_linear_fit!(graph, exp, discr, y_cref;
                                droptol=simlist[j].opt_kwargs[:droptol],
                                linlsqr=simlist[j].opt_kwargs[:linlsqr]                                )
                yy1=get_coeffs(graph,y_cref);
                println("Linear fit degopt (change: $(Float64(norm(yy0-yy1))) )");
                sim=deepcopy(simlist[j]);

            elseif (x=="q");
                println("Quit!");
                break;
            else
                println("Unknown command $x")
                skip=true;
            end
            if (!skip)
                graphlist[j+1]=graph;
                simlist[j+1]=sim;
                commandlist[j+1]=x;
                j=j+1;
            end

        end


        if (x=="s")
            (graph,cref)=runsim(simlist[j-1],graphlist[j-1],cref);
            #
            graphlist[j]=graph;
            simlist[j]=sim;
        end


    end

    return (graph,simlist[1:j],graphlist[1:j],join(commandlist[1:j]))

end



function scale_and_square_degopt(fname,state,kickstart=1)
    graph_org=import_compgraph(fname);;
    (graph,cref)=gen_degopt_by_squaring(graph_org);
    (x,y)=get_degopt_crefs(graph);
    if (kickstart>0)
        x[end][1][1] = 0.001 # Kickstart
    end
    if (kickstart>1)
        x[1][1][1] = 0.001 # Kickstart
    end
    (graph2,cref)=gen_degopt_poly(x,y);

    graph2=Compgraph(Complex{BigFloat},graph2)

    y_cref=get_degopt_crefs(count(values(graph.operations) .== :mult))[2];

    discr=get_disc_discr(state);
    #discr=state.rho*exp.(1im*range(0,2*pi,length=sim.n)[1:end-1]);
    discr=convert.(Complex{BigFloat},discr);
    opt_linear_fit!(graph2, exp, discr, y_cref;
                    errtype = :relerr,
                    linlsqr = :real_svd,
                    droptol = 1e-17)
    graph=deepcopy(graph2);
    return (graph,cref)
end
