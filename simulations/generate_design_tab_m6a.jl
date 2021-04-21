include("reset_all.jl");
include("simulationtools.jl");
include("add_squaring.jl");



m=6;

rho=2.22;  # SID rho m=6
target=Simulation(m,n=200,f=exp,rho=rho)

its=5;
droptol0=1e-10;

base_sim=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
         init=:taylor,
         opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                 :linlsqr => :real_svd,:maxit => its))




sim0_mono21=deepcopy(base_sim);
sim0_mono21.init=:prev;
sim0_mono21.graph=:prev;
sim0_mono21.opt_kwargs[:droptol]=1e-8;
sim0_mono21.opt_kwargs[:γ0]=0.5;
sim0_mono21.opt_kwargs[:maxit]=3;

mono21_init=deepcopy(base_sim);
mono21_init.init=:taylor;
mono21_init.graph=:mono2;
mono21_init.opt_kwargs[:droptol]=1e-10;
mono21_init.opt_kwargs[:γ0]=0.5;
mono21_init.opt_kwargs[:maxit]=0;
(graph_mono21,simlist,graphlist,commandlist)=
        interactive_simulations(mono21_init,sim0_mono21,"IssssssDddddsdddsssssssssssssdddddssssssdddsssssdsGGGssddssssddssssdsddddddgggsddDDDDsssssssssssssssssddDDDDssddssq");
