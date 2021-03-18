# Boostingthecomputationofthematrixexponential

export gen_dic_exp;

"""
    (graph,cref)=gen_dic_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference. The coefficients are directly copied from the paper.

Reference:

* Boosting the computation of the matrix exponential, Sastre, Ibanez, Defez, Appl. Math. Computation, 340, 2019, 206-220.

    """
function gen_dic_exp(k;T=Float64)

    if (k==4)

        c1=convert(T,4.018761610201036*1e-4)
        c2=convert(T,2.945531440279683*1e-3)
        c3=convert(T,-8.709066576837676*1e-3)
        c4=convert(T,4.017568440673568*1e-1)
        c5=convert(T,3.230762888122312*1e-2)
        c6=convert(T,5.768988513026145*1e0)
        c7=convert(T,2.338576034271299*1e-2)
        c8=convert(T,2.381070373870987*1e-1)
        c9=convert(T,2.224209172496374*1e0)
        c10=convert(T,-5.792361707073261*1e0)
        c11=convert(T,-4.130276365929783*1e-2)
        c12=convert(T,1.040801735231354*1e1)
        c13=convert(T,-6.331712455883370*1e1)
        c14=convert(T,3.484665863364574*1e-1)
        c15=convert(T,1)
        c16=convert(T,1)


        # B2=A^2
        v1a=[0;1.0];
        v1b=[0;1.0];

        # B3=y02=A^2*(c1*A^2+c2*A)
        v2a=[0;0;1];
        v2b=[0;c2;c1];

        # B4=(c4*A+c3*A^2+y02)*(c5*A^2+y02)   # First term in y12
        #       Note: y12=B4+c6*y02+c7*A^2
        v3a=[0;c4;c3;1.0];
        v3b=[0;0;c5;1.0];


        # B5=(c9*A+c8*A^2+y12)*(c11*A+c10*y02+y12)    # first term in y22
        #   =(c9*A+(c8+c7)*A^2+c6*y02+B4)*(c11*A+c7*A^2+(c10+c6)*y02+B4)
        v4a=[0;c9;c8+c7;c6;1.0];
        v4b=[0;c11;c7;(c10+c6);1.0];


        # add the additional terms to B5

        y=[c16;c15;(c14+c12*c7);(c13+c12*c6);c12;1.0];

        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_general_poly_recursion(xv,y);


    else
        error("Not implemented for k=$k");
    end


    return (graph,cref);
end
