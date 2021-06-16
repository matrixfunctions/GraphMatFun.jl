# State of the
mutable struct State
    f   # function
    discr
    eltype
    #
    graph
    cref
    #
    params # Dict with keywords
end
function State(f,eltype)
    return State(f,[],eltype,Nothing,Nothing,Dict{Symbol,Any}(:opt_kwargs => Dict()));
end

function scale_and_square!(state)
    degopt=Degopt(state.graph)
    scale!(degopt,1/2)
    square!(degopt)
    (g,c)=gen_degopt_poly(degopt)
    state.graph=g
    state.cref=c
end


function get_disc_discr(state;n=state.params[:n])
    rho=state.params[:rho]
    nn=round(Int,n);
    discr=rho*exp.(1im*range(0,2*pi,length=nn)[1:end-1]);
    discr=convert.(state.eltype,discr);
    return discr;
end
function get_taylor(eltype,n)
    c=1 ./factorial.(convert.(big(real(eltype)),0:(n-1)));
    c=convert.(eltype,c);
end


# Corresponds to a call to opt_gaussnewton!
struct OptSimulation
    kwargs
    cref
end

function OptSimulation(s::State)
    return OptSimulation(s.params[:opt_kwargs],s.cref);
end

# Levenberg-Marquard simulations using LsqFit.curve_fit
struct LMSimulation
    kwargs
end

struct NoOp
end

function LMSimulation(s::State)
    kwargs=Dict();
    if (haskey(s.params,:curve_fit_kwargs))
        kwargs=s.params[:curve_fit_kwargs]
    end
    return LMSimulation(kwargs)
end


# Corresponds to changing a parameter
struct ModifyParam
    param
    factor
end

# Corresponds to a change of all variables
struct KickIt
    factor::Number
    mode
end
function KickIt(state::State,factor)

    if haskey(state.params,:kickit_mode)
        mode=state.params[:kickit_mode]
    else
        mode=0;
    end
    @show mode
    return KickIt(factor,mode);

end


# Initialize State with a graph based on number of multiplications.
function init_state_mult(f,rho,n,mult; eltype=Complex{BigFloat})
    init_state_mult(nothing,f,rho,n,mult; eltype=eltype);
end

function init_state_mult(graphname,f,rho,n,mult; eltype=Complex{BigFloat})
    if (f != exp)
        error("Only exp implemented");
    end
    state=State(f,eltype);
    state.params[:n]=n;
    state.params[:target_n]=n;
    state.params[:rho]=rho;
    state.params[:target_n]=n;

    state.discr=get_disc_discr(state);


    init_state_mult!(state,graphname,mult);
    return state;
end

function init_state_mult!(state,graphname,mult;showmeta=false)

    if (!isnothing(graphname))
        print("Generating $graphname ");
        if (graphname == :sid)
            (graph,cref)=gen_sid_exp(mult)
        elseif (graphname == :bbc)
            (graph,cref)=gen_bbc_basic_exp(m)
        elseif (graphname == :sastre)
            method=:y1s;
            if (m==6)
                method=:h2m
            elseif (m==8)
                method=:z1ps
            end
            (graph,cref)=gen_sastre_basic_exp(m,method)
        end

        if (graphname in [:ps,:mono])

            taylorcoeffs=get_taylor(state.eltype,200);

            # Automatically detect degree
            graph_mult = 0;
            deg0 = m-1;
            graph=nothing; oldgraph=nothing;
            cref=nothing; oldcref=nothing;
            while graph_mult<mult+1
                deg0 += 1;
                c=taylorcoeffs[1:(deg0+1)];
                oldgraph=graph;
                oldcref=cref;
                if (graphname == :ps)
                    (graph,cref)=gen_ps_degopt(c)
                elseif (graphname == :mono)
                    (graph,cref)=gen_monomial_degopt(c)
                end
                graph_mult=count(values(graph.operations) .==:mult);
            end
            # We go one step too far in the automatic determination so
            # use the old graph
            graph=oldgraph;
            cref=oldcref;

            println(" degree $(deg0-1)");
        end

        state.params[:graphname]=string(graphname);
        graph = Compgraph(state.eltype,graph);
        state.graph=graph;
        state.cref=cref;
        if (showmeta)
            showerr(state,showmeta=true)
        end

    end
    return state;
end
function init_state_file!(state,graphname,filename;showmeta=true,scale_and_square=false);

    state.graph=import_compgraph(filename);
    state.params[:graphname]=string(graphname);

    if (scale_and_square)
        scale_and_square!(state)
    else
        (_,cref)=gen_degopt_poly(Degopt(state.graph))
        state.cref=cref;
    end


    if (showmeta)
        println("Generating $graphname based on $filename");
        showerr(state,showmeta=true)
    end

    return state;
end

# initcommand(control::String,state)  --> command object
function initcommand(control,state)
    if (control == "s")
        return OptSimulation(state);
    elseif (control == "l")
        return LMSimulation(state);
    elseif (control == "g")
        return ModifyParam(:γ0,0.8);
    elseif (control == "G")
        return ModifyParam(:γ0,1/0.8);
    elseif (control=="n")
        return ModifyParam(:n,1/2);
    elseif (control=="N")
        return ModifyParam(:n,2);
    elseif (control=="D")
        return ModifyParam(:droptol,2);
    elseif (control=="d")
        return ModifyParam(:droptol,1/2);
    elseif (control=="k")
        return KickIt(state,1e-9);
    elseif (control=="K")
        return KickIt(state,1e-6);
    else
        return NoOp();
    end
end

# runcommand(commandobj, state) -> State
function runcommand(s::OptSimulation,state)
    # Run Gauss-Newton with params in s.params[:opt_kwargs]
    f=state.f;
    graph=deepcopy(state.graph);
    discr=state.discr;
    cref=s.cref;

    opt_gauss_newton!(graph,f,discr;logger=1,
                      stoptol=1e-30,cref=cref,
                      s.kwargs...);

    state=deepcopy(state);

    state.graph=graph # Set the output
    currentname=state.params[:graphname]
    if (!contains(currentname,"+GN"))
        state.params[:graphname]=currentname*"+GN";
    end

    return state;
end

function runcommand(s::ModifyParam,state)
    # Update a parameter some opt_kwargs treated separately
    state=deepcopy(state);
    param = s.param;
    factor = s.factor;
    if (param in [:droptol,:γ0])
        o=state.params[:opt_kwargs]
    else
        o=state.params;
    end
    T=typeof(o[param])
    o[param]=o[param]*factor;
    println("$param=$(o[param])");
    if (param == :n)
        state.discr=get_disc_discr(state);
    end

    return state;
end

function curve_fit_model(x,coeffs,graph,cref)
    p=size(x,1)/2;
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    z=eval_graph(graph,x);
    #set_coeffs!(graph,org_coeffs,cref);
    return [real.(z);imag.(z)];
end


function curve_fit_jacobian_model(x,coeffs,graph,cref)
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    J=eval_jac(graph,x,cref);
    #set_coeffs!(graph,org_coeffs,cref);
    return [real.(J);imag.(J)]
end

function runcommand(s::LMSimulation,state)

    discr=state.discr;
    T=real(state.eltype);
    xdata=complex(T).(discr);
    ydata_c=(state.f.(discr));
    ydata=T.([real.(ydata_c);imag(ydata_c)])

    graph=Compgraph(T,state.graph);
    cref=state.cref;
    p0=get_coeffs(graph,cref);


    model = (x,coeffs) -> curve_fit_model(x,coeffs,graph,cref)
    jacobian_model =
        (x,coeffs) -> curve_fit_jacobian_model(x,coeffs,graph,cref)

    fit = curve_fit(model, jacobian_model, xdata, ydata, p0,#show_trace=true,
                    x_tol=big(1e-20), g_tol=1e-20, maxIter=10,
                    s.kwargs...)

    state.graph=graph # Set the output
    @show typeof(graph)
    currentname=state.params[:graphname]
    if (!contains(currentname,"+LM"))
        state.params[:graphname]=currentname*"+LM";
    end
    return state;

end
function runcommand(s::NoOp,state)
    println("No op");
    return state;
end

function runcommand(s::KickIt,state)
    factor=s.factor;
    (x,y)=get_degopt_crefs(state.graph);
    Random.seed!(0);
    if (s.mode == 0)
        cref=[x[end][1][1]];
        vals=get_coeffs(state.graph,cref);
        vals=vals.*(1 .+ randn(size(vals)));
        set_coeffs!(state.graph,vals,cref);

    elseif (s.mode == 1)
        cref=[x[end][1][1]];
        vals=get_coeffs(state.graph,cref);
        vals=vals .+ randn(size(vals));
        set_coeffs!(state.graph,vals,cref);
    elseif (s.mode == 2)

        cref=[x[1][2][1]];
        vals=get_coeffs(state.graph,cref);
        vals .+= 0.0001
        set_coeffs!(state.graph,vals,cref);
    else

        println("Unknown Kickit Mode");
    end

    return state;
end



function run_sequence(state,predefsims="");
    max_j=1000;
    statelist=Vector{Any}(undef,max_j);
    commandlist=Vector{String}(undef,max_j);
    statelist[1]=state;
    pre_commandlist=split(predefsims,"");

    err=showerr(state,output=false);
    println("Initialerr $(state.params[:graphname]):$err");
    j=1;
    command="";
    while (command != "q")
        state=statelist[j] # Base it on the old state
        print("Commandlist: "*String(join(commandlist[1:(j-1)]))* " ");

        x="";
        if (j <= size(pre_commandlist,1)) # It's stored in
            x=pre_commandlist[j];
            pre_commandlist[j] = "";
        end
        if x == ""
            x=String(read_one_key());
        end

        commandlist[j]=x;
        println("$x");
        if (x=="b") # Back
            j=j-2;
        elseif (x == "q")
            break
        else
            command = initcommand(x,state)
            newstate = runcommand(command,state)
            statelist[j+1] = newstate;
        end

        showerr(statelist[j+1]);

        j += 1;

    end

    state=statelist[j];
    return (state,statelist[1:j],join(commandlist[1:j]))

end


function showmetainfo(state)
    println("rho: $(state.params[:graphname]): ρ=$(state.params[:rho]) ");
end

function showerr(state;output=true,showmeta=false,n=state.params[:target_n])
    f=state.f;
    discr=get_disc_discr(state,n=n);


    graph=state.graph;

    err=norm((f.(discr)-eval_graph(graph,discr))./f.(discr),Inf)
    imagnorm=Float64(norm(imag.(get_coeffs(graph))));
    if (output)
        if (showmeta)
            showmetainfo(state)
        end

        if (imagnorm>0)
            imagstr="norm(imag(coeffs))= $imagnorm"
        else
            imagstr="";
        end
        println("Target error: $(Float64(err)) $imagstr");
    end

    return err
end


function read_one_key(; io = stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io.handle, raw)
    setraw!(true)
    s=read(io, 1)
    setraw!(false)
    return s
end
