
include("simulationtools.jl");

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





m=4;

rho=6.9e-1;  # SID rho m=4
target=Simulation(m,n=200,f=exp,rho=rho)

sid=Simulation(m,n=50,f=exp,rho=rho,init=:lsqr,graph=:sid,
                opt_kwargs=Dict(:logger=> 1))
(sid_graph,_,_)=initsim(sid,0,0);
showerr(target,sid_graph)



droptol0=1e-10;
its=5;



sim0=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))

sim_agressive=deepcopy(sim0);
sim_agressive.opt_kwargs[:γ0]=1.0;
sim_agressive.opt_kwargs[:droptol]=droptol0/100;

sim_agressive2=deepcopy(sim0);
sim_agressive2.opt_kwargs[:γ0]=1.0;
sim_agressive2.opt_kwargs[:droptol]=droptol0/1000;




sim_very_agressive=make_agressive(make_agressive(make_accurate(sim_agressive2,target)));

sim_very_agressive2=make_agressive(sim_very_agressive);



init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:ps,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],10)...,
         "Agressive", repeat([sim_agressive],16)...,
         "Agressive 2", repeat([sim_agressive2],16)...
         ]

(graph1,cref)=run_simlist(simlist)


init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:lsqr,graph=:ps,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],10)...,
         "Agressive", repeat([sim_agressive],16)...,
         "Agressive 2", repeat([sim_agressive2],16)...
         ]
(graph2,cref)=run_simlist(simlist)



init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:mono,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],5)...,
         "Accurate", repeat([make_accurate(sim0,target)],3)...]
(graph3,cref)=run_simlist(simlist);



init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:lsqr,graph=:mono,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],5)...,
         "Accurate", repeat([make_accurate(sim0,target)],3)...]
(graph4,cref)=run_simlist(simlist);


init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:horner,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],30)...,
         "Accurate", repeat([make_accurate(sim0,target)],3)...]
(graph5,cref)=run_simlist(simlist);



init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:horner,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))
simlist=[init_sim,repeat([sim0],30)...,
         "Accurate", repeat([make_accurate(sim0,target)],3)...]
(graph6,cref)=run_simlist(simlist);


init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:bbc,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))

simlist=[init_sim,repeat([sim0],1)...,
         "Agressive", repeat([make_agressive(sim0)],5)...,
         "Accurate", repeat([make_agressive(make_agressive(make_agressive(sim0)))],15)...]
(graph7,cref)=run_simlist(simlist);

err1=showerr(target,graph1,false)
err2=showerr(target,graph2,false)
err3=showerr(target,graph3,false)
err4=showerr(target,graph4,false)
err5=showerr(target,graph5,false)
err6=showerr(target,graph6,false)
err7=showerr(target,graph7,false)

err8=showerr(target,sid_graph,false)
bbc_graph=gen_bbc_basic_exp(4)[1]
err9=showerr(target,bbc_graph,false)
println("PS taylor $err1")
println("PS lsqr $err2")
println("mono taylor $err3")
println("mono lsqr $err4")
println("horner taylor $err5")
println("horner lsqr $err6")
println("BBC + gn $err7")
println("SID $err8")
println("BBC $err9")
