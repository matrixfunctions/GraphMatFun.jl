
using GenericSVD,LinearAlgebra

function get_taylor(n)
    c=1 ./factorial.(0:(n-1));
end
function get_lsqr(n,discr)
    Z=zeros(big(eltype(discr)),size(discr,1),n);
    Z[:,1]=ones(size(discr,1));
    for s=2:n
        Z[:,s]=(Z[:,s-1]).*(discr);
    end
    c=Z\exp.(big.(discr));
    return convert.(eltype(discr),c);
end

m=5;
deg0=m+1;
n=100;
rho=1.0;
discr=rho*exp.(1im*range(0,2*pi,length=n)[1:end-1]);
bdiscr=big.(discr);
c=get_lsqr(deg0+1,bdiscr)
#c=get_taylor(deg0+1)
(graph,cref)=gen_horner_recursive(c);

if (m != count(values(graph.operations) .==:mult))
    @show m, count(values(graph.operations) .==:mult)
    error("Incorrect number of multiplications")
end


# First optimize in ComplexF64
graph=Compgraph(ComplexF64,graph);
opt_gauss_newton!(graph,exp,discr,logger=1,
                  stoptol=1e-16,cref=cref,γ0=0.05,errtype=:relerr);

# First optimize in Complex{BigFloat}
graph=Compgraph(Complex{BigFloat},graph);
# svd and droptol important
opt_gauss_newton!(graph,exp,bdiscr,logger=1,linlsqr=:svd,droptol=1e-10,
                  stoptol=1e-20,cref=cref,γ0=0.5,errtype=:relerr);
@show norm(exp.(bdiscr)-eval_graph(graph,bdiscr),Inf)
