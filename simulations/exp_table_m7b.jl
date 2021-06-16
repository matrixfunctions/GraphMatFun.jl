using GraphMatFun, LsqFit, LinearAlgebra, GenericSVD,Random
include("new_reset_all.jl");
#include("simulationtools.jl");
include("newsimulationtools.jl");


f=exp
m=7
target_n=200 # Measure error wrt this discr
n=50;        # Optimize wrt this discr
rho=6.0;  # Increased SID rho m=7


its=5;
droptol0=1e-10;
opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                :linlsqr => :real_svd,:maxit => its)
base = init_state_mult(f,rho,n,m; eltype=Complex{BigFloat})
base.params[:n]=n;
base.params[:target_n]=target_n;
base.params[:opt_kwargs]=opt_kwargs;
base.params[:kickit_mode]=1;


oldrho="3_59";
oldrho2="2_7";
# First the optimization free
filename="simulations/newgraphs/exp_m$(m)_sid_$oldrho.cgr";
sid_org=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)

#filename="simulations/newgraphs/exp_m$(m)_bbc_$oldrho.cgr";
#bbc_org=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)

filename="simulations/newgraphs/exp_m$(m)_sastre_$oldrho.cgr";
sastre_org=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true)


#sastre_org.cref=get_degopt_crefs(sastre_org.graph);


filename="simulations/newgraphs/exp_m$(m)_ps_$oldrho.cgr";
ps_org=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)


filename="simulations/newgraphs/exp_m$(m)_mono_$oldrho.cgr";
mono_org=init_state_file!(deepcopy(base),:mono,filename,showmeta=true)


filename="simulations/newgraphs/exp_m$(m-1)_mono+GN_$(oldrho2).cgr";
mono0=init_state_file!(deepcopy(base),:mono,filename,showmeta=true,scale_and_square=true)
(mono,simlist,commandlist)=
     run_sequence(mono0,"ddddddddddddddddddddddddddsssGGGsssdddsddsssssdddsssdddsssssq");

#filename="simulations/newgraphs/exp_m$(m)_ps+GN_$oldrho.cgr";
#ps0=init_state_file!(deepcopy(base),:ps,filename,showmeta=true)
#(ps,simlist,commandlist)=
#     run_sequence(ps0,
#     "sq");

filename="simulations/newgraphs/exp_m$(m)_sastre+GN_$oldrho.cgr";
sastre0=init_state_file!(deepcopy(base),:sastre,filename,showmeta=true)
(sastre,simlist,commandlist)=
     run_sequence(sastre0,
                  "ddddddddddsdddddddddddssGGGssddddsdddddsdsssssdddssssddssssddsgggsssssGGGssq");

filename="simulations/newgraphs/exp_m$(m)_sid+GN_$oldrho.cgr";
sid0=init_state_file!(deepcopy(base),:sid,filename,showmeta=true)
(sid,simlist,commandlist)=
     run_sequence(sid_org,
     "")

#filename="simulations/newgraphs/exp_m$(m)_bbc+GN_$oldrho.cgr";
#bbc0=init_state_file!(deepcopy(base),:bbc,filename,showmeta=true)
#(bbc,simlist,commandlist)=
#     run_sequence(bbc0,
#     "sdddddsdddddssssssddddsssssssddddsssssssssssssdddssssssssssddddssssssssdsggggggsdssssssssssssssssGGGGGsssGGGsssq")
#

include("new_print_all.jl");
println("Run include(save_all.jl) if you want to save");
