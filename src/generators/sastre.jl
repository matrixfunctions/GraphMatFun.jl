
export gen_sastre_basic_exp, gen_sastre_basic

"""
    (graph,cref)=gen_sastre_basic_exp(k)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference.

Reference:

*  Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010
    """
function gen_sastre_basic_exp(k)
    if (k==3)
        b=1 ./factorial.(0:8)
        return gen_sastre_basic(b)
    elseif (k==4)
        e0 = 5.018851944498568
        e = [e0, NaN, 0.03806343180936604, 0.017732587443103232]
        c = [0.002193172316532563, 0.0002741465395665704, 4.569108992776174e-5]
        d = [1.3093238729699403, 0.1955094199013519, 0.01626158346315151]
        f = [1.0, 1.0, 0.5, 0.1168293067115003]

        s=3
        p=s
        return gen_sastre_degopt(s,p,c,d,e,f,NaN)
    elseif (k==8) # Table 7
        c10=-6.140022498994532e-17
        c9=-9.210033748491798e-16
        c8=-1.980157255925737e-14
        c7=-4.508311519886735e-13
        c6=-1.023660713518307e-11
        d5=-1.227011356117036e-10
        d4=-6.770221628797445e-9
        d3=-1.502070379373464e-7
        d2=-3.013961104055248e-6
        d1=-5.893435534477677e-5
        e5=-3.294026127901678e-10
        e4=-2.785084196756015e-9
        e3=-4.032817333361947e-8
        e2=-5.100472475630675e-7
        e0=-1.023463999572971e-3
        f5=4.024189993755686e-13
        f4=7.556768134694921e-12
        f3=1.305311326377090e-10
        f2=2.087675698786810e-9
        f1=2.505210838544172e-8
        f0=2.755731922398589e-7

        c=[c6;c7;c8;c9;c10]
        d=[d1;d2;d3;d4;d5]
        e=[e0;NaN;e2;e3;e4;e5]

        f=[f0;f1;f2;f3;f4;f5]

        a=1 ./factorial.(0:9)

        s=5
        p=10
        return gen_sastre_degopt(s,p,c,d,e,f,a)
    else
        error("Not implemented k=$k")
    end
end


"""
    (graph,cref)=gen_sastre_basic(b)

Computes the degree-8 polynomial

    p(z)=b[1]+z*b[2]+z^2*b[3]+...+z^8*b[9]

according to Example 3.1 in the reference.

Reference:

*  Efficient evaluation of matrix polynomials, J. Sastre. Linear Algebra and its Applications ,Volume 539, 2018, Pages 229-250, https://doi.org/10.1016/j.laa.2017.11.010
    """
function gen_sastre_basic(b)
# Equations (16) - (32)
    if (size(b,1) !=9)
        error("Not implemented");
    end

    b0=b[1]
    b1=b[2]
    b2=b[3]
    b3=b[4]
    b4=b[5]
    b5=b[6]
    b6=b[7]
    b7=b[8]
    b8=b[9]


    c4=sqrt(b8) # plus minus?
    c3=b7/(2*c4)
    d2_plus_e2=(b6-c3^2)/c4
    d1=(b5-c3*d2_plus_e2)/c4


    e2_num_sqrt=(d1-(c3/c4)*d2_plus_e2)^2+4*(c3/c4)*(b3+(c3^2/c4)*d1-(c3/c4)*b4)

    e2_num=(c3/c4)*d2_plus_e2-d1+sqrt(e2_num_sqrt)  #plus minus
    e2=e2_num/(2*c3/c4)
    d2=d2_plus_e2-e2
    f2=b2
    f1=b1
    f0=b0


    e0=(b3-d1*e2)/c3 # Not explicitly documented?




    e=[e0;NaN;e2]
    c=[c3;c4]
    f=[f0;f1;f2]
    d=[d1;d2]

    s=2
    p=s
    return gen_sastre_degopt(s,p,c,d,e,f,NaN)
    #gen_degopt_poly(x,z);


end


# Internal use only
# Transforms formula (34)-(35) + (52) to degopt form
# Input are the coefficients given in the paper
# s int
# d=[d1,...ds]
# c=[c_(s+1)...c_(2*s)]  # size = s
# e=[e0;NaN;e2;...e_s]  # size = s+1
# f=[f0;...f_s] # size = s+1
# a=[a0;...a_(p-1)] # size = p   Can be NaN if p=s
# Nof mult: s+1+p/s = s+1+v
function gen_sastre_degopt(s,p,c,d,e,f,a)
    T=eltype(c)
    x = Vector{Tuple{Vector{T},Vector{T}}}()

    for j=1:s-1
        push!(x,(vcat(zero(T),one(T),zeros(T,j-1))
                ,vcat(zeros(T,j),one(T)))
             )
    end

    # y0s
    push!(x,(vcat(zeros(T,s),one(T))
            ,vcat(zero(T),c[1:s]))
         )

    # first term y1s
    push!(x,(vcat(zero(T),d[1:s],one(T))
            ,vcat(zero(T),zero(T),e[3:(s+1)],one(T)))
         )

    # y1s
    ys1 = vcat(f[1:s+1],e[1],one(T))
    if (p==s) # Only (32)-(34)

        z=ys1
        return gen_degopt_poly(x,z)
    else #Apply PS-scheme as in (52)
        v = convert(Int,p/s) # Assumed to be an integer, according to paper, p = v*s

        # Evaluate y1s*x^s
        push!(x,(ys1
               ,vcat(zeros(T,s),one(T),zeros(T,2)))
             )

         # Evaluate rest of the multiplications with PS-scheme
         for i = reverse(0:v-2)
             lower_idx = (i+1)*s+1
             upper_idx = lower_idx+s-1
             idx = lower_idx:upper_idx
             c = view(a, idx)
             xl = vcat(c,zero(T),zeros(T,v-i),one(T)) # Coeffs for P_{v-i} + A^s*P_{v+1-i}
             xr = vcat(zeros(T,s),one(T),zeros(T,v+1-i)) # Coeffs for A^s
             push!(x, (xl,xr)) # A^s*( P_{v-i} + A^s*P_{v+1-i} )
         end

        z = vcat(a[1:s],zeros(T,2+v),one(T))
        return gen_degopt_poly(x,z)
    end



end
