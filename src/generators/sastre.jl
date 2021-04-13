
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

    elseif (k==4)
        e0 = 5.018851944498568
        e = [e0, NaN, 0.03806343180936604, 0.017732587443103232]
        c = [0.002193172316532563, 0.0002741465395665704, 4.569108992776174e-5]
        d = [1.3093238729699403, 0.1955094199013519, 0.01626158346315151]
        f = [1.0, 1.0, 0.5, 0.1168293067115003]

        s=3;
        return gen_sastre_degopt(s,c,d,e,f)
    elseif (k==6)
        c10=-6.140022498994532E-17
        c9=-9.210033748491798E-16
        c8=-1.980157255925737E-14
        c7=-4.508311519886735E-13
        c6=-1.023660713518307E-11
        d5=-1.227011356117036E-10
        d4=-6.770221628797445E-9
        d3=-1.502070379373464E-7
        d2=-3.013961104055248E-6
        d1=-5.893435534477677E-5
        e5=-3.294026127901678E-10
        e4=-2.785084196756015E-9
        e3=-4.032817333361947E-8
        e2=-5.100472475630675E-7
        e0=-1.023463999572971E-3
        f5=4.024189993755686E-13
        f4=7.556768134694921E-12
        f3=1.305311326377090E-10
        f2=2.087675698786810E-9
        f1=2.505210838544172E-8
        f0=2.755731922398589E-7

        c=[c6;c7;c8;c9;c10];
        d=[d1;d2;d3;d4;d5];
        e=[e0;NaN;e2;e3;e4;e5];

        f=[f0;f1;f2;f3;f4;f5];

        s=5;
        return gen_sastre_degopt(s,c,d,e,f)
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




    e=[e0;NaN;e2];
    c=[c3;c4];
    f=[f0;f1;f2];
    d=[d1;d2];

    s=2;
    return gen_sastre_degopt(s,c,d,e,f)
    #gen_degopt_poly(x,z);


end


# Internal use only
# Transforms formula (34)-(35) to degopt form
# Input are the coefficients given in the paper
# s int
# d=[d1,...ds]
# c=[c_(s+1)...c_(2*s)]  # size = s
# e=[e0;NaN;e2;...e_s]  # size = s+1
# f=[f0;...f_s] # size = s+1
# Nof mult: s+1
function gen_sastre_degopt(s,c,d,e,f)
    T=eltype(c)
    x = Vector{Tuple{Vector{T},Vector{T}}}()

    for j=1:s-1
        push!(x,([0.0;1.0;zeros(T,j-1)],[zeros(T,j);1.0]));
    end

    # y0s
    push!(x,([zeros(T,s);1.0],[0.0;c[1:s]]));

    # first term y1s
    push!(x,([0.0;d[1:s];1.0],[0.0;0.0;e[3:(s+1)];1.0]));

    # y1s
    z=[f[1:s+1];e[1];1.0];

    @show x
    @show z
    return gen_degopt_poly(x,z);
end
