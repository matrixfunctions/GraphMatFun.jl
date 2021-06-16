include("reset_all.jl");
include("simulationtools.jl");




m=5;

rho=1.9;  # SID rho m=5
target=Simulation(m,n=200,f=exp,rho=rho)

its=5;
droptol0=1e-10;

base_sim=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
         init=:taylor,
         opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                 :linlsqr => :real_svd,:maxit => its))


degopt_cref=vec_degopt(get_degopt_coeffs(m));


sid=deepcopy(base_sim);
sid.graph=:sid;
(sid_org,_,_)=initsim(sid,0,0);
showerr(target,sid_org)


#sastre=deepcopy(base_sim);
#sastre.graph=:sastre;
#(sastre_org,_,_)=initsim(sastre,0,0);
#showerr(target,sastre_org)

bbc=deepcopy(base_sim);
bbc.graph=:bbc;
(bbc_org,_,_)=initsim(bbc,0,0);
showerr(target,bbc_org)



sim0=deepcopy(base_sim);
sim0.init=:prev;
sim0.graph=:prev;



ps1_init=deepcopy(base_sim);
ps1_init.init=:taylor;
graph_org=import_compgraph("simulations/graphs/exp_m5_PS_taylor_1_68.cgr");
ps1_init.graph=(graph_org,degopt_cref)
(graph_ps1,simlist,graphlist,commandlist)=
interactive_simulations(ps1_init,sim0,
                        "IsnddddssddsddsddslssssddsddsddslsssddsddsddsNsddsssssddsddssslssssssssdsq");
# Err: 1.2E-11





mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
mono1_init.graph=:mono;
graph_org=import_compgraph("simulations/graphs/exp_m5_mono_taylor_1_68.cgr");
mono1_init.graph=(graph_org,degopt_cref)
(graph_mono1,simlist,graphlist,commandlist)=
        interactive_simulations(mono1_init,sim0,"IssddsddsddsddsddssNslssslsdgggggsssssssssssssssssssssGGsGsssssssssssssddssssssssssssssssssdsddsddssssssssssssssssssdsddsssssssssssssssssssddsssssGsssssssssslsssssq");
# Err: 1.1E-11



sastre_init=deepcopy(base_sim);
sastre_init.init=:taylor;
graph_org=import_compgraph(
    "simulations/graphs/exp_m5_Sastre+_1_68.cgr")
sastre_init.graph=(graph_org,degopt_cref);
(graph_sastre,simlist,graphlist,commandlist)=
        interactive_simulations(sastre_init,sim0,"Isddsddsddsddsddsddsddsddsddsssssq");
#

bbc_init=deepcopy(base_sim);
bbc_init.init=:taylor;
bbc_init.graph=:bbc;
(graph_bbc,simlist,graphlist,commandlist)=
        interactive_simulations(bbc_init,sim0,"Isssssddssssddsssssssssssdddssdddssssdddssdddssssssssssdggdrq");
#



sid_init=deepcopy(base_sim);
sid_init.init=:taylor;
graph_org=import_compgraph(
    "simulations/graphs/exp_m5_SID+_1_68.cgr")
sid_init.graph=(graph,degopt_cref);
sid_init.graph=:sid;
(graph_sid,simlist,graphlist,commandlist)=
        interactive_simulations(sid_init,sim0,"IsssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddsssssssssdddssssddsdssssssrq");


include("print_all.jl");
println("Run include(save_all.jl) if you want to save");
