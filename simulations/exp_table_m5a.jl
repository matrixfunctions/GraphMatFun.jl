using GraphMatFun, LsqFit, LinearAlgebra
include("new_reset_all.jl");
#include("simulationtools.jl");
include("newsimulationtools.jl");


f=exp
m=5
target_n=200 # Measure error wrt this discr
n=50;        # Optimize wrt this discr
rho=1.68;  # SID rho m=5


its=5;
droptol0=1e-10;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=n_target;
base.params[:opt_kwargs]=opt_kwargs;


# First the optimization free
sid_org=init_state_mult!(deepcopy(base),:sid,m,showmeta=true)

bbc_org=init_state_mult!(deepcopy(base),:bbc,m,showmeta=true)

filename="simulations/newgraphs/exp_m4_sastre+GN_0_69.cgr"
sastre_org=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true,
                            scale_and_square=true)


#sastre_org.cref=get_degopt_crefs(sastre_org.graph);


ps_org=init_state_mult!(deepcopy(base),:ps,m,showmeta=true)


mono_org=init_state_mult!(deepcopy(base),:mono,m,showmeta=true)



(mono,simlist,commandlist)=
     run_sequence(mono_org,"snsssssssssddssddsssssssssssssdsddGGssssq");

(ps,simlist,commandlist)=
     run_sequence(ps_org,
     "snsssssssssddssddsssssssssssssdssddssssssssssssssddsddssNssNssdssdsssddssssssssssssssssssssddssssssddsssssssdddsssssssNssdssssssssddssssssssssdssssq");



(sastre,simlist,commandlist)=
     run_sequence(sastre_org,
                  "sssnssNdddddssdddssGGssssdddssssddsssssddsddssddsssssq");


(sid,simlist,commandlist)=
     run_sequence(sid_org,
     "sssssddssssddsssssssssssdddssdddssssdddssdddsssNsssdddsssssssssdddssssddsdssssssq")

(bbc,simlist,commandlist)=
     run_sequence(bbc_org,
     "sssssddssssddsssssssssssdddssdddssssdddssdddssssssssssdggdq")


include("new_print_all.jl");
println("Run include(save_all.jl) if you want to save");
