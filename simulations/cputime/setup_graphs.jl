using MAT,LinearAlgebra;

graph_m6=import_compgraph("simulations/graphs/exp_m6_mono_taylor_2_7.cgr");;
graph_m6_SID=import_compgraph("simulations/graphs/exp_m6_SID_2_22.cgr");;
graph_m7=import_compgraph("simulations/graphs/exp_m7_SID+_3_59.cgr");;

graph_m6=Compgraph(Float64,graph_m6);
graph_m6_SID=Compgraph(Float64,graph_m6_SID);
compress_graph!(graph_m6_SID);

graph_m7=Compgraph(Float64,graph_m7);


#A=matread("simulations/cputime/n2000_2_2.mat")["A"];
n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*I;
A0=2.5*A0/8
# opnorm(A0,1)=1;


(graph_native,_)=gen_exp_native_jl(A0);

names=["m6_mono_taylor_2_7","m6_SID_2_22","m7_SIDplus_3_59","native_jl"];
graphs=[graph_m6,graph_m6_SID,graph_m7,graph_native]


priohelp=Dict();
addit=-10;
priohelp[:Ba2]=-10+addit;
priohelp[:Bb2]=-10+addit;
priohelp[:Ba3]=-9+addit;
priohelp[:Ba3_2]=-9+addit;
priohelp[:Bb3]=-8+addit;
priohelp[:Bb3_2]=-8+addit;

priohelp[:Ba4]=-7+addit;
priohelp[:Ba4_2]=-7+addit;
priohelp[:Ba4_3]=-7+addit;
priohelp[:Bb4]=-6+addit;
priohelp[:Bb4_2]=-6+addit;
priohelp[:Bb4_3]=-6+addit;


priohelp[:Ba5]=-5+addit;
priohelp[:Ba5_2]=-5
priohelp[:Ba5_3]=-5+addit;
priohelp[:Ba5_4]=-5+addit;
priohelp[:Bb5]=-4+addit;
priohelp[:Bb5_2]=-4+addit;
priohelp[:Bb5_3]=-4+addit;
priohelp[:Bb5_4]=-4+addit;



priohelp[:Ba6]=-3+addit;
priohelp[:Ba6_2]=-3
priohelp[:Ba6_3]=-3+addit;
priohelp[:Ba6_4]=-3+addit;
priohelp[:Ba6_5]=-3+addit;
priohelp[:Bb6]=-2+addit;
priohelp[:Bb6_2]=-2+addit;
priohelp[:Bb6_3]=-2+addit;
priohelp[:Bb6_4]=-2+addit;
priohelp[:Bb6_5]=-2+addit;


priohelp[:Ba7]=-1.5+addit;
priohelp[:Ba7_2]=-1.5
priohelp[:Ba7_3]=-1.5+addit;
priohelp[:Ba7_4]=-1.5+addit;
priohelp[:Ba7_5]=-1.5+addit;
priohelp[:Ba7_6]=-1.5+addit;
priohelp[:Bb7]=-1.1+addit;
priohelp[:Bb7_2]=-1.1+addit;
priohelp[:Bb7_3]=-1.1+addit;
priohelp[:Bb7_4]=-1.1+addit;
priohelp[:Bb7_5]=-1.1+addit;
priohelp[:Bb7_6]=-1.1+addit;
