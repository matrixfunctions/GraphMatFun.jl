
# Code related to the paper
#
# Efficient evaluation of matrix polynomials
# J. Sastre. Linear Algebra and its Applications
# Volume 539, 15 February 2018, Pages 229-250
# https://doi.org/10.1016/j.laa.2017.11.010


"""
    (graph,cref)=gen_sastre_basic_exp(k)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference.

Reference:

*  Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010
    """
function gen_sastre_basic_exp(k)
    graph=Compgraph(Float64);
    cref=[];
    if (k==3)
        b=1 ./factorial.(0:8);
        return gen_sastre_basic(b)
    else
        error("Not implemented k=$k")
    end




end


function gen_sastre_basic(b)

    @show b
    if (size(b,1) !=9)
        error("Not implemented");
    end

    b0=b[1];
    b1=b[2];
    b2=b[3];
    b3=b[4];
    b4=b[5];
    b5=b[6];
    b6=b[7];
    b7=b[8];
    b8=b[9];


    c4=sqrt(b8); # plus minus?
    c3=b7/(2*c4)
    d2_plus_e2=(b6-c3^2)/c4;
    d1=(b5-c3*d2_plus_e2)/c4;


    e2_num_sqrt=(d1-(c3/c4)*d2_plus_e2)^2+4*(c3/c4)*(b3+(c3^2/c4)*d1-(c3/c4)*b4);

    e2_num=(c3/c4)*d2_plus_e2-d1+sqrt(e2_num_sqrt);  #plus minus
    e2=e2_num/(2*c3/c4);
    d2=d2_plus_e2-e2;
    f2=b2;
    f1=b1;
    f0=b0;


    e0=(b3-d1*e2)/c3; # Not explicitly documented?


    # Now build by transforming (14) to degopt form

    T=Float64;
    x = Vector{Tuple{Vector{T},Vector{T}}}(undef,3)

    x[1]=([0, 1], [0, 1.0]);
    x[2]=([0, 0, 1], [0, c3, c4]) # y02
    x[3]=([0, d1, d2, 1], [0, 0, e2, 1]); # first term in y12
    z=[f0, f1, f2, e0, 1];


    return gen_degopt_poly(x,z);


end
