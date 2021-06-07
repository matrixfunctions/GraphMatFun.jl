include("reset_all.jl");
include("simulationtools.jl");
include("add_squaring.jl");



m=8;

rho=6.5*2;  # 2* theta for m =7
rho=13.5;
target=Simulation(m,n=200,f=exp,rho=rho)

its=5;
droptol0=1e-10;

base_sim=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
         init=:taylor,
         opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                 :linlsqr => :real_svd,:maxit => its))


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


println("Mono21");
mono21_init=deepcopy(base_sim);
mono21_init.init=:taylor;
graph_org=import_compgraph("simulations/graphs/exp_m7_mono_taylor_6_0.cgr");;
(graph,cref)=gen_degopt_by_squaring(graph_org,kickstart=2);
(graph,cref)=gen_degopt_poly(Degopt(graph));
mono21_init.graph=(graph,cref);
mono21_init.opt_kwargs[:droptol]=1e-20;
mono21_init.opt_kwargs[:γ0]=0.5;
mono21_init.opt_kwargs[:maxit]=0;
mono21_sim0=deepcopy(sim0);
mono21_sim0.opt_kwargs[:droptol]=1e-15

(graph_mono1,simlist,graphlist,commandlist)=
interactive_simulations(mono21_init,mono21_sim0,"IsssssdddsssssddddssssssddddddsssssssdddsssssssssssdddssssGGsssddddggggggssssssssGGGsssGGGssssdddddggggggggsssssssGGGsGGGGGssssdddddggggggggssssssssGGGsGGGGGsssssddddgggggggsssGGGGGGGssssssdddggggggggggssssssggssssssGsGsGsGsGGsGGsGGGsGGsss");

                        "IsssssdddsssssddddssssssddddddsssssssdddsssssssssssdddssssGGsssddddggggggssssssssGGGsssGGGssssdddddggggggggsssssssGGGsGGGGGssssdddddggggggggssssssssGGGsGGGGGsssssddddgggggggsssGGGGGGGssssssdddddggggggggsssssssq")

"IsssssdddsssssddddssssssddddddsssssssdddsssssssssssdddssssGGsssddddggggggssssssssGGGsssGGGssssdddddggggggggsssssssGGGsGGGGGssssdddddggggggggssssssssGGGsGGGGGss"

#                        "IsssssdddsssssddddssssssddddddsssssssdddssssssssssddsssssssssdddssrrddsggsssGGsssssssssssddsdggssssssGGssssssssdgggsgggsdssssssssssssssssGGsssssGGGssssssssssdddgggssssssggsssssssGGssggssssssssssssssssssGsGssssssGsGsGsGssssssssNssssGssssq");
asd
ps1_init=deepcopy(base_sim);
ps1_init.init=:taylor;
ps1_init.graph=:ps;
(graph_ps1,simlist,graphlist,commandlist)=
interactive_simulations(ps1_init,sim0,
                        "IssnssrsddsddsddddsdddddsrsssdddssrsssssddddsssssssdddsdsdggsssNsssNssdddsssssssdddsssssssssssssssssssdddsssssssddddssssssddddssddddsssssGsGsssssddsssssssssrssrssrsrq");
# Err: 1.2E-11

ps2_init=deepcopy(base_sim);
ps2_init.init=:lsqr;
ps2_init.graph=:ps;
(graph_ps2,simlist,graphlist,commandlist)=
interactive_simulations(ps2_init,sim0,
                        "IsnNsnssrsddsddsddddsdddddsrsssdddssrsssssddddsssssssdddsdsdggsssNsssNssdddsssssssdddsssssssssssssssssssdddsssssssddddssssssddddssddddsssssGsGsssssddsssssssssrssrssrsNssssssrq");

# Err: 6E-12



mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
mono1_init.graph=:mono;
(graph_mono1,simlist,graphlist,commandlist)=
        interactive_simulations(mono1_init,sim0,"IsnsddddssdddssdddsssdddsssdddsssdddssssssssssrsssddsddssssssssssssdddsssssssssssssddssrsddsssssssssssssssNsssssssssssdsddssssssssssssrsddsddsddsddsssssrssssssssssddsrssddsddsddsddsssdddsddsrq");
# Err: 1.1E-11

mono2_init=deepcopy(base_sim);
mono2_init.init=:lsqr;
mono2_init.graph=:mono;
(graph_mono2,simlist,graphlist,commandlist)=
        interactive_simulations(mono2_init,sim0,"IsnsddddssdddssdddsssdddsssdddsssdddssssssssssrsssddsddssssssssssssdddsssssssssssssddssrsddsssssssssssssssNsssssssssssdsddssssssssssssrsddsddsddsddsssssrq");
# Err: 1.1e-11



#sastre_init=deepcopy(base_sim);
#sastre_init.init=:taylor;
#sastre_init.graph=:sastre;
#(graph_sastre,simlist,graphlist,commandlist)=
#        interactive_simulations(sastre_init,sim0,"Isssssddssssddssssssssssssddsssssssddsssddsssssssssssssssssssssssssdssssssssssssssssddsrssssssrssssssdsdssssrsssssssdssssssssssdsssrssssrq");
##  1.448620852496371e-15
#

#bbc_init=deepcopy(base_sim);
#bbc_init.init=:taylor;
#bbc_init.graph=:bbc;
#(graph_bbc,simlist,graphlist,commandlist)=
#        interactive_simulations(bbc_init,sim0,"Isssssddssssddsssssssssssdddssdddssssdddssdddsssrq");
## 2.2665827420785147e-14
#
#

sid_init=deepcopy(base_sim);
sid_init.init=:taylor;
sid_init.graph=:sid;
(graph_sid,simlist,graphlist,commandlist)=
        interactive_simulations(sid_init,sim0,"IsnddddddddddddddssssssssssNsssddsdddddddsddddddddddddddsddssdddddddsssssssssddddddddsssssssssssssddssssssssssssssrdddsssssssssssssddsssssssssssssrq");



include("print_all.jl");
println("Run include(save_all.jl) if you want to save");
