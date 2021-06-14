include("simulationtools.jl");
using LsqFit
using Random;

m=4;

rho=6.9e-1;  # SID rho m=4
target=Simulation(m,n=100,f=exp,rho=rho)

discr=rho*exp.(2im*big(pi)*(1:target.n)/target.n)
its=5;
droptol0=1e-10;

base_sim=Simulation(m,n=50,f=exp,rho=rho,eltype=Complex{BigFloat},
         init=:taylor,
         opt_kwargs=Dict(:logger=> 0,:γ0 => 0.5,:droptol => droptol0,
                 :linlsqr => :real_svd,:maxit => its))
sim0=deepcopy(base_sim);


#(graph_ps1,cref)=gen_ps_degopt(1 ./ factorial.(Int128.(0:9)));
(graph_ps1,cref)=gen_ps_degopt(1 ./ factorial.(Int128.(0:12)));
(graph_ps1,_)=gen_degopt_poly(normalize!(Degopt(graph_ps1)))
#cref=cref[5:end];
cref=cref[20:end];
#cref=cref[11:end];

count( values(graph_ps1.operations) .== :mult)
#cref=cref[[1:10;end]];
graph=deepcopy(graph_ps1);
showerr(sim0,graph)

T=BigFloat
graph=Compgraph(T,graph)


function model(x,coeffs)
    p=size(x,1)/2;
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    z=eval_graph(graph,x);
    set_coeffs!(graph,org_coeffs,cref);
    return [real.(z);imag.(z)];
end


function jacobian_model(x,coeffs)
    org_coeffs=get_coeffs(graph,cref);
    set_coeffs!(graph,coeffs,cref)
    J=eval_jac(graph,x,cref);
    set_coeffs!(graph,org_coeffs,cref);
    return [real.(J);imag.(J)]
end

xdata=complex(T).(discr);
ydata_c=(exp.(discr));
ydata=T.([real.(ydata_c);imag(ydata_c)])
Random.seed!(0);
p0=get_coeffs(graph,cref).*(1 .+ 0.02*randn(size(cref,1)));
#p0=get_coeffs(graph,cref);

set_coeffs!(graph,p0,cref);
showerr(sim0,graph)


fit = curve_fit(model, jacobian_model, xdata, ydata, p0,show_trace=true,
                x_tol=big(1e-20), g_tol=1e-20, maxIter=50000)

set_coeffs!(graph,fit.param,cref);
showerr(sim0,graph)

#p0=get_coeffs(graph,cref);
#
#fit = curve_fit(model, jacobian_model, xdata, ydata, p0)
#
#showerr(sim0,graph)
