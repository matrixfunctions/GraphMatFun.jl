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


#monosqr1_init=deepcopy(base_sim);
#monosqr1_init.init=:taylor;
#graph_org=import_compgraph("simulations/graphs/exp_m6_mono_taylor_2_7.cgr");;
#showerr(target,graph_org)
#(graph,cref)=gen_degopt_by_squaring(graph_org);
#(x,y)=get_degopt_coeffs(graph);
#x[end][1][1] = 0.001 # Kickstart
##y[1] += 1e-6;
#(graph,cref)=gen_degopt_poly(x,y);
#graph=Compgraph(Complex{BigFloat},graph);
#y_cref=get_degopt_coeffs(count(values(graph.operations) .== :mult))[2];
#showerr(target,graph)
#discr=sim0.rho*exp.(1im*range(0,2*pi,length=sim0.n)[1:end-1]);
#discr=convert.(Complex{BigFloat},discr);
#
#opt_linear_fit!(graph, exp, discr, y_cref;
#                  errtype = :relerr,
#                  linlsqr = :real_backslash,
#                  droptol = 1e-17)
#showerr(target,graph)
#
#
##(xx,yy)=get_degopt_coeffs(graph);
#
#monosqr1_init.graph=(graph,cref);
#monosqr1_init.opt_kwargs[:droptol]=1e-14;
#monosqr1_init.opt_kwargs[:γ0]=0.5;
#monosqr1_init.opt_kwargs[:maxit]=0;
#(graph_monosqr1,simlist,graphlist,commandlist)=
#  interactive_simulations(monosqr1_init,sim0,"IssssdddddsssGGssssssrsssrssssssrssssssdddsssrssssssdddssssddddssssssdddddddssrssssssssddddddddssssss");
##IssssdddddsssGGssssssrsssrssssssrssssssdddsssrssssssdddssssddddssssssdddddddssrssssssssddddddddssssssddddsssssssdsrq");
### rho 6 => 2.6422405729310504e-14
##
## IssssdddddsssGGssssssrsssrssssssrssssssdddsssrssssssdddssssddddssssss => 2.3E-14
## IssssdddddsssGGssssssrsssrssssssrssssssdddsssrssssssdddssssddddssssssdddddddssrssssssssddddddddssssss => 7.6E-16
##monosqr1_init.graph=(graph_monosqr1,cref);
##monosqr1_init.rho=5.5; sim0.rho=5.5;
##(graph_monosqr1b,simlist,graphlist,commandlist)=
##        interactive_simulations(monosqr1_init,sim0)
##
##asd

println("PS1");
ps1_init=deepcopy(base_sim);
ps1_init.init=:taylor;
ps1_init.graph=:ps;
(graph,cref)=scale_and_square_degopt(
    "simulations/graphs/exp_m6_PS_taylor_2_7.cgr",sim0,0)
ps1_init.graph=(graph,cref);
(graph_ps1,simlist,graphlist,commandlist)=
interactive_simulations(ps1_init,sim0,
                        "IsssdddddddddddsssddddssddddssdddddddsssddssssdlsssssssssssssssssssssssGgggggssssssssssGsssssssssssssssssssssssssssssssssddsq");
#                        "IsnsssssssdssGGGsssssdddsssdddNsssssssssssdsssssssddssssdsssddsssddssddsssrsssssssdssdsrsdddsdggssssssssssddssddssssssssssddssssssssddssddssgdsssrssssssssrrrdggssssssssssssrssssssssssssssGsssssddssddssrq");
#                        "IssnssrsddsddsddddsdddddsrsssdddssrsssssddddsssssssdddsdsdggsssNsssNssdddsssssssdddsssssssssssssssssssdddsssssssddddssssssddddssddddsssssGsGsssssddsssssssssrssrssrsrq");
# Err: 1.2E-11

#ps2_init=deepcopy(base_sim);
#ps2_init.init=:lsqr;
#ps2_init.graph=:ps;
#(graph_ps2,simlist,graphlist,commandlist)=
#interactive_simulations(ps2_init,sim0,
#                        "IsnNsnssrsddsddsddddsdddddsrsssdddssrsssssddddsssssssdddsdsdggsssNsssNssdddsssssssdddsssssssssssssssssssdddsssssssddddssssssddddssddddsssssGsGsssssddsssssssssrssrssrsNssssssrq");
#
## Err: 6E-12
#



#(graph_mono1,simlist,graphlist,commandlist)=
#        interactive_simulations(mono1_init,sim0,"IsnsssssssdssGGGsssssdddsssdddNsssssssssssdsssssssddssssdsssddsssddssddsssrsssssssdssdsrsdddsdggssssssssssddssddssssssssssddssssssssddssddssgdsssrssssssssrrrdggssssssssssssrssssssssssssssGsssssddssddssrq");
#
println("MONO1");
mono1_init=deepcopy(base_sim);
mono1_init.init=:taylor;
(graph,cref)=scale_and_square_degopt(
    "simulations/graphs/exp_m6_mono_taylor_2_7.cgr",sim0,1)
mono1_init.graph=(graph,cref);
mono1_init.opt_kwargs[:droptol]=1e-14;
mono1_init.opt_kwargs[:γ0]=0.5;
mono1_init.opt_kwargs[:maxit]=0;
mono1_sim0=deepcopy(sim0);
#mono1_sim0.opt_kwargs[:linlsqr]=:svd;
(graph_mono1,simlist,graphlist,commandlist)=
interactive_simulations(mono1_init,sim0,
                        "IssssdddddsssGGssssssrsssrssssssrssssssdddsssrssssssdddssssddddssssssdddddddssrssssssssddddddddssssssq");






#mono2_init=deepcopy(base_sim);
#mono2_init.init=:lsqr;
#mono2_init.graph=:mono;
#(graph_mono2,simlist,graphlist,commandlist)=
#        interactive_simulations(mono2_init,sim0,"IsnsssssssdssGGGsssssdddsssdddNsssssssssssdsssssssddssssdsssddsssddssddsssrsssssssdssdsrsdddsdggssssssssssddssddssssssssssddssssssssddssddssgdsssrssssssssrrrdggssssssssssssrssssssssssssssGsssssddssddssrq");
#
#


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



println("SID+: Using previous");
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

#sid_init=deepcopy(base_sim);
#sid_init.init=:taylor;
##sid_init.graph=:sid;
#(graph,cref)=scale_and_square_degopt(
#    "simulations/graphs/exp_m6_SID+_2_7.cgr",sim0,1)
#sid_init.graph=(graph,cref);
#(graph_sid,simlist,graphlist,commandlist)=
#        interactive_simulations(sid_init,sim0,"IsnddddddddddddddssssssssssNsssddsdddddddsddddddddddddddsddssdddddddsssssssssddddddddsssssssssssssddssssssssssssssrdddsssssssssssssddsssssssssssssrq");
#
#

include("print_all.jl");
println("Run include(save_all.jl) if you want to save");
