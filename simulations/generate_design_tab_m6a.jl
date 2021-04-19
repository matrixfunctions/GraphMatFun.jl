include("reset_all.jl");
include("simulationtools.jl");




m=6;

rho=2.22;  # SID rho m=6
target=Simulation(m,n=200,f=exp,rho=rho)

its=5;
droptol0=1e-10;

base_sim=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
         init=:taylor,
         opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                 :linlsqr => :svd,:maxit => its))


sid=deepcopy(base_sim);
sid.graph=:sid;
(sid_org,_,_)=initsim(sid,0,0);
showerr(target,sid_org)

#sastre=deepcopy(base_sim);
#sastre.graph=:sastre;
#(sastre_org,_,_)=initsim(sastre,0,0);
#showerr(target,sastre_org)
#



sim0=deepcopy(base_sim);
sim0.init=:prev;
sim0.graph=:prev;


ps1_init=deepcopy(base_sim);
ps1_init.init=:taylor;
ps1_init.graph=:ps;
(graph_ps1,simlist,graphlist,commandlist)=
interactive_simulations(ps1_init,sim0,
                        "IsnsssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddssNssNssdsrsdsssddssssssssrssssssssssssddssssssddssssssrsdddsssssssNssdssssssssddssssssssrssdssssssssssrdssggsssssssssssrssssssssssssrssssrq");
# Err: 1.2E-11

ps2_init=deepcopy(base_sim);
ps2_init.init=:lsqr;
ps2_init.graph=:ps;
(graph_ps2,simlist,graphlist,commandlist)=
interactive_simulations(ps2_init,sim0,
                        "IsnsssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddssNssNssdsrsdsssddssssssssrssssssssssssddssssssddssssssrsdddsssssssNssdssssssssddssssssssrssdssssssssssrdssggsssssssssssrssssssssssssrssssrq");

# Err: 6E-12



mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
mono1_init.graph=:mono;
(graph_mono1,simlist,graphlist,commandlist)=
        interactive_simulations(mono1_init,sim0,"IsnsssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddssNssNssdsrsdsssddssssssssrssssssssssssddssssssddssssssrsdddsssssssNssdssssssssddssssssssrssdssssssssssrdssggsssssssssssrssssssssssssrssssrq");
# Err: 1.1E-11

mono2_init=deepcopy(base_sim);
mono2_init.init=:lsqr;
mono2_init.graph=:mono;
(graph_mono2,simlist,graphlist,commandlist)=
        interactive_simulations(mono2_init,sim0,"IsnsssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddssNssNssdsrsdsssddssssssssrssssssssssssddssssssddssssssrsdddsssssssNssdssssssssddssssssssrssdssssssssssrdssggsssssssssssrssssssssssssrssssrq");
# Err: 1.1e-11


sid_init=deepcopy(base_sim);
sid_init.init=:taylor;
sid_init.graph=:sid;
(graph_sid,simlist,graphlist,commandlist)=
        interactive_simulations(sid_init,sim0,"IsssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddsssssssssdddssssddsdssssssrq");



include("print_all.jl");
println("Run include(save_all.jl) if you want to save");
