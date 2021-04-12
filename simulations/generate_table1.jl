
include("simulationtools.jl");



m=4;

rho=6.9e-1;  # SID rho m=4
target=Simulation(m,n=200,f=exp,rho=rho)

sid=Simulation(m,n=50,f=exp,rho=rho,init=:lsqr,graph=:sid,
                opt_kwargs=Dict(:logger=> 1))
(sid_graph,_,_)=initsim(sid,0,0);
showerr(target,sid_graph)



droptol0=1e-10;
its=5;


init_sim=Simulation(m,n=100,f=exp,rho=rho,init=:taylor,graph=:ps,
                   eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 1.0,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))
#                opt_kwargs=Dict(:logger=> 1,:maxit => 40))

sim0=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
                opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                                :linlsqr => :svd,:maxit => its))

sim_agressive=deepcopy(sim0);
sim_agressive.opt_kwargs[:γ0]=1.0;
sim_agressive.opt_kwargs[:droptol]=droptol0/100;

sim_agressive2=deepcopy(sim0);
sim_agressive2.opt_kwargs[:γ0]=1.0;
sim_agressive2.opt_kwargs[:droptol]=droptol0/1000;


sim_smooth=deepcopy(sim0);
sim_smooth.opt_kwargs[:γ0]=0.5;
sim_smooth.opt_kwargs[:droptol]=droptol0;

sim_smooth2=deepcopy(sim0);
sim_smooth2.opt_kwargs[:γ0]=0.5;
sim_smooth2.opt_kwargs[:droptol]=droptol0/100;


sim_very_agressive=make_agressive(make_agressive(make_accurate(sim_agressive2,target)));

sim_very_agressive2=make_agressive(sim_very_agressive);



simlist=Vector();
simlist=[init_sim,sim0,sim0,sim0,sim0,sim0,sim0,sim0,sim0,sim0,sim0,
         "Agressive", sim_agressive,sim_agressive,sim_agressive,sim_agressive,sim_agressive,sim_agressive,sim_agressive,sim_agressive,
         "Very agressive", sim_agressive2, sim_agressive2, sim_agressive2, sim_agressive2, sim_agressive2, sim_agressive2, sim_agressive2,
         "Final", sim_very_agressive, sim_very_agressive, sim_very_agressive, sim_very_agressive, sim_very_agressive, sim_very_agressive,
         "Final", sim_very_agressive2, sim_very_agressive2]

#push!(simlist,init_sim);
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,sim0)
#push!(simlist,"Agressive")
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,sim_agressive)
#push!(simlist,"Very agressive")
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,sim_agressive2)
#push!(simlist,make_accurate(sim_agressive2,target))
#push!(simlist,"Final");
#push!(simlist,sim_very_agressive)
#push!(simlist,sim_very_agressive)
#push!(simlist,sim_very_agressive)
#push!(simlist,sim_very_agressive)
#push!(simlist,sim_very_agressive2)
#push!(simlist,sim_very_agressive2)
##
#push!(simlist,"maxit:$(sim.opt_kwargs[:maxit])");
#push!(simlist,"droptol:$(sim.opt_kwargs[:droptol])");
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#sim=deepcopy(sim);
#sim.n=200
#sim.opt_kwargs[:maxit]=5
#push!(simlist,"n:$(sim.n)");
#push!(simlist,"maxit:$(sim.opt_kwargs[:maxit])");
#push!(simlist,"droptol:$(sim.opt_kwargs[:droptol])");
#push!(simlist,sim);
#sim=deepcopy(sim);
#sim.rho=rho
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#sim=deepcopy(sim);
#sim.opt_kwargs[:droptol]=1e-14
#push!(simlist,"droptol:$(sim.opt_kwargs[:droptol])");
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);
#sim=deepcopy(sim);
#sim.opt_kwargs[:droptol]=1e-16
#sim.opt_kwargs[:maxit]=32
#push!(simlist,"maxit:$(sim.opt_kwargs[:maxit])");
#push!(simlist,"droptol:$(sim.opt_kwargs[:droptol])");
#push!(simlist,sim);
#push!(simlist,sim);
#push!(simlist,sim);

#
#sim=deepcopy(sim);
#sim.rho=0.69;
#sim.opt_kwargs[:maxit]=10;
#push!(simlist,sim);
#
#sim=deepcopy(sim);
#sim.rho=0.69;
#push!(simlist,sim);
#
#
#sim=deepcopy(sim);
#push!(simlist,sim);
#
#sim=deepcopy(sim)
#sim.opt_kwargs[:maxit]=20;
#sim.opt_kwargs[:γ0]=0.3;
#push!(simlist,sim);
#push!(simlist,sim);
#
#graph=runsim(sim1,0,0);

graph=0; cref=[];
for sim=simlist
    global graph,cref
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
