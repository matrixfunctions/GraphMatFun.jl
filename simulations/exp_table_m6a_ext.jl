using GraphMatFun, LsqFit, LinearAlgebra, GenericSVD,Random
include("new_reset_all.jl");
#include("simulationtools.jl");
include("newsimulationtools.jl");


f=exp
m=6
target_n=200 # Measure error wrt this discr
n=50;        # Optimize wrt this discr
rho=2.7;  # SID


its=5;
droptol0=1e-10;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=target_n;
base.params[:opt_kwargs]=opt_kwargs;
base.params[:kickit_mode]=1;

# First the optimization free
sid_org=init_state_mult!(deepcopy(base),:sid,m,showmeta=true)

#bbc_org=init_state_mult!(deepcopy(base),:bbc,m,showmeta=true)

#filename="simulations/newgraphs/exp_m4_sastre+GN_0_69.cgr"
sastre_org=init_state_mult!(deepcopy(base),:sastre,m,showmeta=true)


#sastre_org.cref=get_degopt_crefs(sastre_org.graph);


ps_org=init_state_mult!(deepcopy(base),:ps,m,showmeta=true)


mono_org=init_state_mult!(deepcopy(base),:mono,m,showmeta=true)


(mono,simlist,commandlist)=
     run_sequence(mono_org,"sknsssssssssddssddsssssssssssssdsddGGssssdddsGssdddsssdddssssddsssssssssssssddsssssddddsddddssssq");

(ps,simlist,commandlist)=
     run_sequence(ps_org, "kssnsKsssdddddsdddsGGsssGssdddssssssddssdddddsssddddsssddddsssssssssssddddddssGgsssssssq");


(sastre,simlist,commandlist)=
     run_sequence(sastre_org,
                  "sssnssNdddddddddddddssGGGsssddddddddsssdgggggssssssGGssssssssdsssGGssssddsddsddsssssssssq");


(sid,simlist,commandlist)=
     run_sequence(sid_org,
     "sssnssNdddddddddddddssGGGsssddsddsddsddsddddsdddsddddddssssddsssssddsssssddssssssq")

#(bbc,simlist,commandlist)=
#     run_sequence(bbc_org,
#     "sssssddssssddsssssssssssdddssdddssssdddssdddssssssssssdggdq")
#

include("new_print_all.jl");
println("Run include(save_all.jl) if you want to save");
