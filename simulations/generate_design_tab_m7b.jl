include("reset_all.jl");
include("simulationtools.jl");



include("add_squaring.jl");

m=7;

rho=6.0;  # SID rho m=7
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
#



sim0=deepcopy(base_sim);
sim0.init=:prev;
sim0.graph=:prev;


ps1_init=deepcopy(base_sim);
ps1_init.init=:taylor;
graph_org=import_compgraph("simulations/graphs/exp_m7_PS_taylor_3_59.cgr");
ps1_init.graph=(graph_org,degopt_cref)
(graph_ps1,simlist,graphlist,commandlist)=
interactive_simulations(ps1_init,sim0,
                        "IsNddddddddddssddsssdddsdddsq");
#                        "IssnssrsddsddsddddsdddddsrsssdddssrsssssddddsssssssdddsdsdggsssNsssNssdddsssssssdddsssssssssssssssssssdddsssssssddddssssssddddssddddsssssGsGsssssddsssssssssrssrssrsrq");
# Err: 1.2E-11

ps2_init=deepcopy(base_sim);
ps2_init.init=:lsqr;
graph_org=import_compgraph("simulations/graphs/exp_m7_PS_lsqr_3_59.cgr");
ps2_init.graph=(graph_org,degopt_cref)
(graph_ps2,simlist,graphlist,commandlist)=
interactive_simulations(ps2_init,sim0,
                        "IsNddddddddddssddsssdddsdddsdddddssdddsssdddddsslssddddsdddddslssssddddddddssslsssssss");

# Err: 6E-12



mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
mono1_init.graph=:mono;
graph_org=import_compgraph("simulations/graphs/exp_m7_mono_taylor_3_59.cgr");
(x,y)=get_degopt_coeffs(graph_org);
x[end][1][1] = 0.001 # Kickstart
x[end-1][2][1] = 0.001 # Kickstart
#y[1] += 1e-6;
(graph_org,cref)=gen_degopt_poly(x,y);

mono1_init.graph=(graph_org,degopt_cref)
(graph_mono1,simlist,graphlist,commandlist)=
        interactive_simulations(mono1_init,sim0,"IslsslssssslddddddddssssssssddddsslslsssssdddddddssssssssddssdssslssssssssssddsssddsslssGssssddddssDdddslsssddsq");





mono2_init=deepcopy(base_sim);
mono2_init.init=:lsqr;
mono2_init.graph=:mono;
graph_org=import_compgraph("simulations/graphs/exp_m7_mono_lsqr_3_59.cgr");
(x,y)=get_degopt_coeffs(graph_org);
x[end][1][1] = 0.001 # Kickstart
x[end-1][2][1] = 0.001 # Kickstart
#y[1] += 1e-6;
(graph_org,cref)=gen_degopt_poly(x,y);
mono2_init.graph=(graph_org,degopt_cref)

(graph_mono2,simlist,graphlist,commandlist)=
        interactive_simulations(mono2_init,sim0,"IslsslssssslddddddddssssssssddddsslslsssssdddddddssssssssddssdssslssssssssssddsssddsslssGssssddddssDdddslsssddsq");




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
graph_org=import_compgraph("simulations/graphs/exp_m7_SID+_3_59.cgr");
(x,y)=get_degopt_coeffs(graph_org);
x[end][1][1] = 0.001 # Kickstart
x[end-1][2][1] = 0.001 # Kickstart
#y[1] += 1e-6;
(graph_org,cref)=gen_degopt_poly(x,y);
sid_init.graph=(graph_org,degopt_cref)

(graph_sid,simlist,graphlist,commandlist)=
        interactive_simulations(sid_init,sim0,"IsssdddddddddddddsslssssddddddsssssdddsssssssssdddsssssssGsdddsssssssdddddssslsssdddsssssssdddddslssssssssdddslsssssssddddslsssssssssdddslslsslsssssssGslsNsdsdgslslsssssssssslsddsslssssssslsq");



include("print_all.jl");
println("Run include(save_all.jl) if you want to save");
