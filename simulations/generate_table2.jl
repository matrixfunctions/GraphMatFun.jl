
include("simulationtools.jl");




m=5;

rho=1.68;  # SID rho m=5
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
                        "IsnsssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddssNssNssdsrsdsssddssssssssrssssssssssssddssssssddssssssrsdddsssssssNssdssssssssddssssssssrssdss");
# Err: 1.2E-11

ps2_init=deepcopy(base_sim);
ps2_init.init=:lsqr;
ps2_init.graph=:ps;
(graph_ps2,simlist,graphlist,commandlist)=
interactive_simulations(ps2_init,sim0,
                        "Issssssssssddssddsrssssssssssssdsrsddssssssssssssssddsdrq");
# Err: 6E-12



mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
mono1_init.graph=:mono;
(graph_mono1,simlist,graphlist,commandlist)=
        interactive_simulations(mono1_init,sim0,"Issssssssssddssddsrssssssssssssdsrsddssssssssssssssddsddrq");
# Err: 1.1E-11

mono2_init=deepcopy(base_sim);
mono2_init.init=:lsqr;
mono2_init.graph=:mono;
(graph_mono2,simlist,graphlist,commandlist)=
        interactive_simulations(mono2_init,sim0,"Issssssssssddssddsrssssssssssssdsrsddssssssssssssssddsq");
# Err: 1.1e-11



#sastre_init=deepcopy(base_sim);
#sastre_init.init=:taylor;
#sastre_init.graph=:sastre;
#(graph_sastre,simlist,graphlist,commandlist)=
#        interactive_simulations(sastre_init,sim0,"Isssssddssssddssssssssssssddsssssssddsssddsssssssssssssssssssssssssdssssssssssssssssddsrssssssrssssssdsdssssrsssssssdssssssssssdsssrssssrq");
##  1.448620852496371e-15
#

bbc_init=deepcopy(base_sim);
bbc_init.init=:taylor;
bbc_init.graph=:bbc;
(graph_bbc,simlist,graphlist,commandlist)=
        interactive_simulations(bbc_init,sim0,"Isssssddssssddsssssssssssdddssdddssssdddssdddsssrq");
# 2.2665827420785147e-14



sid_init=deepcopy(base_sim);
sid_init.init=:taylor;
sid_init.graph=:sid;
(graph_sid,simlist,graphlist,commandlist)=
        interactive_simulations(sid_init,sim0,"IsssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddsssssssssdddssssddsdssssssrq");

err1=showerr(target,graph_ps1,false)
err2=showerr(target,graph_ps2,false)
err3=showerr(target,graph_mono1,false)
err4=showerr(target,graph_mono2,false)
#err5=showerr(target,graph_sastre,false)
err6=showerr(target,graph_bbc,false)
err7=showerr(target,graph_sid,false)

#errX1=showerr(target,sastre_org,false)
errX2=showerr(target,sid_org,false)
bbc_graph=gen_bbc_basic_exp(4)[1]
errX3=showerr(target,bbc_graph,false)
println("PS taylor $err1")
println("PS lsqr $err2")
println("mono taylor $err3")
println("mono lsqr $err4")
#println("Sastre+ $err5")
println("BBC+ $err6")
println("SID+ $err7")
#println("Sastre $errX1")
println("SID $errX2")
println("BBC $errX3")
