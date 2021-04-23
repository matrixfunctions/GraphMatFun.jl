using LinearAlgebra;

graph_m6=import_compgraph("simulations/graphs/exp_m6_mono_taylor_2_7.cgr");;
graph_m6_SID=import_compgraph("simulations/graphs/exp_m6_SID_2_22.cgr");;
graph_m7=import_compgraph("simulations/graphs/exp_m7_mono_taylor_6_0.cgr");;
graph_m7_SID=import_compgraph("simulations/graphs/exp_m7_SID+_6_0.cgr");;

#graph_m6=Compgraph(Float64,graph_m6);
#graph_m6_SID=Compgraph(Float64,graph_m6_SID);
#compress_graph!(graph_m6_SID);
#
#graph_m7=Compgraph(Float64,graph_m7);
#

#A=matread("simulations/cputime/n2000_2_2.mat")["A"];
n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*I;
A0[5,3] += 0.0001; # Don't end up in hermitian special case
col=2;
if col==1
    A0=2.5*A0/8
else
    A0=5.5*A0/8
end

# opnorm(A0,1)=1;


(graph_native,_)=gen_exp_native_jl(A0);
(graph_native2,_)=gen_exp_native_jl(ones(1,1)*5.5);

names=["m6_mono_taylor_2_7","m6_SID_2_22","m7_mono_taylor_6_0", "m7_SIDplus_6_0","native_73_jl","native_83_jl"];
graphs=[graph_m6,graph_m6_SID,graph_m7,graph_m7_SID,graph_native,graph_native2]
for (i,g)=enumerate(graphs)
    if (count(values(g.operations) .== :ldiv) == 0)
        degopt=Degopt(g);
        normalize!(degopt);
        (g,_)=gen_degopt_poly(degopt);
    end
    x=get_coeffs(g)
    if (norm(imag.(x))<eps()*100)
        set_coeffs!(g,real.(x));
    else
        error("non-real coeffs");
    end
    g=Compgraph(Float64,g);
    compress_graph!(g);
    graphs[i]=g;
end





priohelp=Dict();
priohelp[:B2]=-1e8;
priohelp[:Ba3_2]=-90;
priohelp[:Ba3]=-90;
priohelp[:Bb3_2]=-90;
priohelp[:Bb3]=-90;
priohelp[:B3]=-90;


priohelp[:Ba4_2]=-80;
priohelp[:Ba4_3]=-80;
priohelp[:Ba4]=-80;
priohelp[:Bb4_2]=-80;
priohelp[:Bb4_3]=-80;
priohelp[:Bb4]=-80;
priohelp[:B4]=-80;

priohelp[:Ba5_2]=-70;
priohelp[:Ba5_3]=-70;
priohelp[:Ba5_4]=-70;
priohelp[:Ba5]=-70;
priohelp[:Bb5_2]=-70;
priohelp[:Bb5_3]=-70;
priohelp[:Bb5_4]=-70;
priohelp[:Bb5]=-70;
priohelp[:B5]=-70;


priohelp[:Ba6_2]=-60;
priohelp[:Ba6_3]=-60;
priohelp[:Ba6_4]=-60;
priohelp[:Bb6_2]=-60;
priohelp[:Bb6_3]=-60;
priohelp[:Bb6_4]=-60;
priohelp[:Ba6]  =-60;
priohelp[:Bb6]  =-60;
priohelp[:B6]  =-60;

priohelp[:Ba7_2]=0;
priohelp[:Ba7_3]=0;
priohelp[:Bb7_2]=0;
priohelp[:Bb7_3]=0;

priohelp[:T2k2]=1000;
priohelp[:T2k3]=1000;
priohelp[:T2k4]=1000;
priohelp[:T2k5]=1000;

#
#priohelp[:Ba4_2]=-90;
#priohelp[:Bb4_2]=-90;
#priohelp[:Ba5_2]=-90;
#priohelp[:Bb5_2]=-90;
#priohelp[:Ba6_2]=-90;
#priohelp[:Bb6_2]=-90;
#priohelp[:Ba7_2]=-90;
#priohelp[:Bb7_2]=-90;
#
#
#priohelp[:Ba3]=-90;
#priohelp[:Bb3]=-90;
#priohelp[:B3]=-90;
#priohelp[:T2k2]=-80;
#priohelp[:T2k3]=-80;
#priohelp[:Ba4_2]=-70;
#priohelp[:Bb4_2]=-70;
#
#priohelp[:Ba4_2]=-60;
#priohelp[:Ba4_3]=-60;
#priohelp[:Ba4]=-60;
#priohelp[:Bb4_2]=-60;
#priohelp[:Bb4_3]=-60;
#priohelp[:Bb4]=-60;
#priohelp[:B4]=-60;
#
#
#priohelp[:Ba5_2]=-50;
#priohelp[:Ba5_3]=-50;
#priohelp[:Bb5_2]=-50;
#priohelp[:Bb5_3]=-50;
#
#priohelp[:Ba6_2]=-50;
#priohelp[:Ba6_3]=-50;
#priohelp[:Bb6_2]=-50;
#priohelp[:Bb6_3]=-50;
#
#
#priohelp[:Ba7_2]=-50;
#priohelp[:Ba7_3]=-50;
#priohelp[:Bb7_2]=-50;
#priohelp[:Bb7_3]=-50;
#
#priohelp[:Ba4]=-60;
#priohelp[:Bb4_2]=-60;
#priohelp[:Bb4_3]=-60;
#priohelp[:Bb4]=-60;
#priohelp[:B4]=-90;
#
#addit=-500;
#factor=10;
#for (i,s)=enumerate(get_topo_order_degopt(6))
#    priohelp[s]=addit+factor*i;
#end
#

#priohelp=Dict();
#addit=-10;
#priohelp[:Ba2]=-10+addit;
#priohelp[:Bb2]=-10+addit;
#priohelp[:B2]=-10+addit;
#priohelp[:Ba3]=-9+addit;
#priohelp[:Ba3_2]=-9+addit;
#priohelp[:Bb3]=-8+addit;
#priohelp[:Bb3_2]=-8+addit;
#priohelp[:B3]=-8+addit;
#
#priohelp[:Ba4]=-7+addit;
#priohelp[:Ba4_2]=-7+addit;
#priohelp[:Ba4_3]=-7+addit;
#priohelp[:Bb4]=-6+addit;
#priohelp[:Bb4_2]=-6+addit;
#priohelp[:Bb4_3]=-6+addit;
#priohelp[:B4]=-6+addit;
#
#priohelp[:Ba5]=-5+addit;
#priohelp[:Ba5_2]=-5
#priohelp[:Ba5_3]=-5+addit;
#priohelp[:Ba5_4]=-5+addit;
#priohelp[:Bb5]=-4+addit;
#priohelp[:Bb5_2]=-4+addit;
#priohelp[:Bb5_3]=-4+addit;
#priohelp[:Bb5_4]=-4+addit;
#priohelp[:B5]=-4+addit;
#
#
#
#priohelp[:Ba6]=-3+addit;
#priohelp[:Ba6_2]=-3
#priohelp[:Ba6_3]=-3+addit;
#priohelp[:Ba6_4]=-3+addit;
#priohelp[:Ba6_5]=-3+addit;
#priohelp[:Bb6]=-2+addit;
#priohelp[:Bb6_2]=-2+addit;
#priohelp[:Bb6_3]=-2+addit;
#priohelp[:Bb6_4]=-2+addit;
#priohelp[:Bb6_5]=-2+addit;
#priohelp[:B6]=-2+addit;
#
#
#priohelp[:Ba7]=-1.5+addit;
#priohelp[:Ba7_2]=-1.5
#priohelp[:Ba7_3]=-1.5+addit;
#priohelp[:Ba7_4]=-1.5+addit;
#priohelp[:Ba7_5]=-1.5+addit;
#priohelp[:Ba7_6]=-1.5+addit;
#priohelp[:Bb7]=-1.1+addit;
#priohelp[:Bb7_2]=-1.1+addit;
#priohelp[:Bb7_3]=-1.1+addit;
#priohelp[:Bb7_4]=-1.1+addit;
#priohelp[:Bb7_5]=-1.1+addit;
#priohelp[:Bb7_6]=-1.1+addit;
#priohelp[:B7]=-1.1+addit;
#
