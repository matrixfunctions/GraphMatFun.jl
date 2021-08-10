export graph_sastre_exp, graph_sastre_poly, graph_sastre_yks_degopt

"""
    (graph,cref)=graph_sastre_exp(k,method=:auto)

Computes a polynomial evaluation approximating the exponential using `k` matrix
multiplications following a `method` given in the reference. The schemes are
embedded into `degop`-format, see [`graph_degopt`](@ref).

The methods are (from the paper referenced below):

  - `:ps_degopt`, Paterson–Stockmeyer method embedded into `degopt`-format.
  - `:y1s`, given by equations (34)-(35)
  - `:z1ps`, given by equations (34)-(35) and (52)
  - `:h2m`, given by equations (34)-(35) and (69)

Not all combinations of `k` and `method` are implemented. Available ones are:

  - `k<3`, `method`=`:ps_degopt`
  - `k=3`, `method`=`:y1s`, as per Table 4 in the reference
  - `k=4`, `method`=`:y1s`
  - `k=6`, `method`=`:h2m`, as per Table 11 in the reference
  - `k=8`, `method`=`:z1ps`, as per Table 7 in the reference

The default option `method=:auto` will choose method according to the value of
`k`, as prescribed above.

**Reference**

1. J. Sastre. "Efficient evaluation of matrix polynomials". Linear Algebra and
   its Applications, 539:229-250, 2018.
   DOI: [10.1016/j.laa.2017.11.010](https://doi.org/10.1016/j.laa.2017.11.010)
"""
function graph_sastre_exp(k, method = :auto)
    if (method == :auto)
        if (k < 3)
            method = :ps_degopt
        elseif (k == 3) || (k == 4)
            method = :y1s
        elseif (k == 6)
            method = :h2m
        elseif (k == 8)
            method = :z1ps
        else
            error("This number of multiplication is not available.")
        end
    end

    if (k < 3) && (method == :ps_degopt)
        if k == 2
            deg = 4
        elseif k == 1
            deg = 2
        else
            deg = 1
        end
        return graph_ps_degopt(1 ./ factorial.(0:deg))

    elseif (k == 3) && (method == :y1s) # Table 4
        e0 = 2.974307204847627
        e2 = 1.225521150112075e-1
        d1 = 8.765009801785554e-1
        d2 = 7.665265321119147e-2
        c3 = 1.992047682223989e-2
        c4 = 4.980119205559973e-3
        f0 = 1.0
        f1 = 1.0
        f2 = 0.5

        c = [c3; c4]
        d = [d1; d2]
        e = [e0; NaN; e2]

        f = [f0; f1; f2]

        s = 2
        p = 0
        return graph_sastre_z1ps_degopt(s, p, c, d, e, f, NaN)

    elseif (k == 4) && (method == :y1s) # Not tabulated?
        e0 = 5.018851944498568
        e = [e0, NaN, 0.03806343180936604, 0.017732587443103232]
        c = [0.002193172316532563, 0.0002741465395665704, 4.569108992776174e-5]
        d = [1.3093238729699403, 0.1955094199013519, 0.01626158346315151]
        f = [1.0, 1.0, 0.5, 0.1168293067115003]

        s = 3
        p = 0
        return graph_sastre_z1ps_degopt(s, p, c, d, e, f, NaN)

    elseif (k == 6) && (method == :h2m) # Table 11 and equation (69)
        e0 = 7.922322450524197
        e2 = 5.317514832355802e-2
        d1 = 2.868706220817633e-1
        d2 = -3.289442879547955e-2
        c3 = 1.052151783051235e-3
        c4 = 4.675683454147702e-4
        f0 = 3.096467971936040
        f1 = 7.723603212944010e-1
        f2 = 1.673139636901279e-1
        c = [c3; c4]
        d = [d1; d2]
        e = [e0; NaN; e2]
        f = [f0; f1; f2]

        e0p = 1.930814505527068
        e2p = 2.771400028062960e-2
        d1p = 3.968985915411500e-1
        d2p = 2.219811707032801e-2
        c3p = 2.688394980266927e-3
        c4p = 4.675683454147702e-4
        f0p = 0.0
        f1p = 2 * 1.614743005681339e-1 #Error in paper? Seems we should use 2*f1'
        f2p = 8.092036376147299e-2
        cp = [c3p; c4p]
        dp = [d1p; d2p]
        ep = [e0p; NaN; e2p]
        fp = [f0p; f1p; f2p]

        s = 2
        graph_sastre_h2m_degopt(s, (c, cp), (d, dp), (e, ep), (f, fp))

    elseif (k == 8) && (method == :z1ps) # Table 7
        c10 = -6.140022498994532e-17
        c9 = -9.210033748491798e-16
        c8 = -1.980157255925737e-14
        c7 = -4.508311519886735e-13
        c6 = -1.023660713518307e-11
        d5 = -1.227011356117036e-10
        d4 = -6.770221628797445e-9
        d3 = -1.502070379373464e-7
        d2 = -3.013961104055248e-6
        d1 = -5.893435534477677e-5
        e5 = -3.294026127901678e-10
        e4 = -2.785084196756015e-9
        e3 = -4.032817333361947e-8
        e2 = -5.100472475630675e-7
        e0 = -1.023463999572971e-3
        f5 = 4.024189993755686e-13
        f4 = 7.556768134694921e-12
        f3 = 1.305311326377090e-10
        f2 = 2.087675698786810e-9
        f1 = 2.505210838544172e-8
        f0 = 2.755731922398589e-7

        c = [c6; c7; c8; c9; c10]
        d = [d1; d2; d3; d4; d5]
        e = [e0; NaN; e2; e3; e4; e5]

        f = [f0; f1; f2; f3; f4; f5]

        a = 1 ./ factorial.(0:9)

        s = 5
        p = 10
        return graph_sastre_z1ps_degopt(s, p, c, d, e, f, a)

    else
        error("Not implemented k=$k and method=$method")
    end
end

"""
    (graph,cref)=graph_sastre_poly(b)

Computes the degree-8 polynomial

    p(z)=b[1]+z*b[2]+z^2*b[3]+...+z^8*b[9]

according to Example 3.1 in the reference.

**Reference**

1. J. Sastre. "Efficient evaluation of matrix polynomials". Linear Algebra and
   its Applications, 539:229-250, 2018. DOI:
   [10.1016/j.laa.2017.11.010](https://doi.org/10.1016/j.laa.2017.11.010)
"""
function graph_sastre_poly(b)
    # Equations (16) - (32)
    if (size(b, 1) != 9)
        error("Not implemented for length(b)=$k")
    end

    b0 = b[1]
    b1 = b[2]
    b2 = b[3]
    b3 = b[4]
    b4 = b[5]
    b5 = b[6]
    b6 = b[7]
    b7 = b[8]
    b8 = b[9]

    c4 = sqrt(b8) # plus minus?
    c3 = b7 / (2 * c4)
    d2_plus_e2 = (b6 - c3^2) / c4
    d1 = (b5 - c3 * d2_plus_e2) / c4

    e2_num_sqrt =
        (d1 - (c3 / c4) * d2_plus_e2)^2 +
        4 * (c3 / c4) * (b3 + (c3^2 / c4) * d1 - (c3 / c4) * b4)

    e2_num = (c3 / c4) * d2_plus_e2 - d1 + sqrt(e2_num_sqrt)  #plus minus
    e2 = e2_num / (2 * c3 / c4)
    d2 = d2_plus_e2 - e2
    f2 = b2
    f1 = b1
    f0 = b0

    e0 = (b3 - d1 * e2) / c3 # Not explicitly documented?

    e = [e0; NaN; e2]
    c = [c3; c4]
    f = [f0; f1; f2]
    d = [d1; d2]

    s = 2
    p = s
    return graph_sastre_z1ps_degopt(s, p, c, d, e, f, NaN)
end

"""
    (graph,cref)=graph_sastre_yks_degopt(k,s,c)

Transforms the polynomial evaluation format given by equations (62)-(65) in the
reference to `degop`-format. The `graph` is a representation of ```y_{k,s}```.
Input `c` is a grouping of the coefficients as given by the representation
(62)-(65). `c` is a `Vector{Vector{Vector}}` of length `k+1`, representing

    [
    [c_i^{(0,1)}, c_i^{(0,2)}]
    [c_i^{(1,1)}, c_i^{(1,2)}, c_i^{(1,3)}, c_i^{(1,4)}, c_i^{(1,5)}, c_i^{(1,6)}]
    ...
    [c_i^{(k,1)}, c_i^{(k,2)}, c_i^{(k,3)}, c_i^{(k,4)}, c_i^{(k,5)}, c_i^{(k,6)}]
    ]

Hence, `c[1]` contains two vectors, the first of length `s`-1 and the second of
length `s`. (Note: In the first vector the constant for I is set to zero) and
`c[2]` up to `c[k+1]` conatins six vectors: Even-numbered vectors have length
`s`+1 Odd-numbered vectors have length j-1, where j is the intex in `c`, e.g.,
`c[2][1]` has one element and `c[3][5]` have two. For example, equations
(57)-(59) are implemented as:

    c = [
    [[c15, c16], [0.0, 0, 0]], # (57)
    [[1.0], [0, c13, c14], [1.0], [c11, 0, c12], [c10], [0.0, 0, 0]], # (58)
    [[0.0, 1], [0, c8, c9], [c7, 1], [0, c6, 0], [c4, c5], [c1, c2, c3]], # (59)
    ]

**Reference**

1. J. Sastre. "Efficient evaluation of matrix polynomials". Linear Algebra and
   its Applications, 539:229-250, 2018. DOI:
   [10.1016/j.laa.2017.11.010](https://doi.org/10.1016/j.laa.2017.11.010)
"""
function graph_sastre_yks_degopt(k, s, c)
    T = eltype(eltype(eltype(c)))
    x = Vector{Tuple{Vector{T},Vector{T}}}()

    for j = 1:s-1
        push!(
            x,
            (vcat(zero(T), one(T), zeros(T, j - 1)), vcat(zeros(T, j), one(T))),
        )
    end

    # y_0s first term
    push!(
        x,
        (vcat(zeros(T, s), one(T)), vcat(zero(T), c[1][1][1:s])), #C_i^{(0,1)}
    )
    if (k == 0)
        z = vcat(c[1][2][1:s+1], one(T)) #C_i^{(0,2)}
        return graph_degopt(x, z)
    else

        # y_js = g_js + sum_i=0^(j-1) c_i^(j,5) y_is + sum_i=0^s c_i^(j,6) x^i
        # Where g_js are the result of the multiplication that is saved. By
        # expanding y_is into similar sums there is a recursion of equivalent
        # adjusted coefficients that needs to be calculated as iterations
        # proceed. Other terms have to adjusted for. Store recursion in vector
        # below.

        # Vector of adjustments of monomial coeffs
        adj_coeffs_mon = [c[1][2][1:s+1]]
        # Vector of adjustments of coeffs for pre-calculated-terms
        adj_coeffs_mult = [zeros(T, 0)]

        # y_1s first term - incorporating second term of y0s
        a = c[2][2][1:s+1] + c[2][1][1] * adj_coeffs_mon[1]
        b = c[2][4][1:s+1] + c[2][3][1] * adj_coeffs_mon[1]
        push!(x, (vcat(a, c[2][1][1]), vcat(b, c[2][3][1])))

        push!(adj_coeffs_mon, c[2][6][1:s+1] + c[2][5][1] * adj_coeffs_mon[1])
        push!(adj_coeffs_mult, [c[2][5][1]])

        for j = 2:k
            # yjs first term - incorporating second and third terms of lower yls
            # for l < j.
            a1 = c[j+1][2][1:s+1]
            b1 = c[j+1][4][1:s+1]
            c1 = c[j+1][6][1:s+1]
            for l = 1:j
                a1 += c[j+1][1][l] * adj_coeffs_mon[l]
                b1 += c[j+1][3][l] * adj_coeffs_mon[l]
                c1 += c[j+1][5][l] * adj_coeffs_mon[l]
            end

            a2 = c[j+1][1][1:j]
            b2 = c[j+1][3][1:j]
            c2 = c[j+1][5][1:j]
            for l = 2:j
                a2[1:l-1] += c[j+1][1][l] * adj_coeffs_mult[l][1:l-1]
                b2[1:l-1] += c[j+1][3][l] * adj_coeffs_mult[l][1:l-1]
                c2[1:l-1] += c[j+1][5][l] * adj_coeffs_mult[l][1:l-1]
            end

            push!(x, (vcat(a1, a2), vcat(b1, b2)))

            push!(adj_coeffs_mon, c1)
            push!(adj_coeffs_mult, c2)
        end

        z = vcat(adj_coeffs_mon[k+1], adj_coeffs_mult[k+1], one(T))
    end

    return graph_degopt(x, z)
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
function graph_sastre_z1ps_degopt(s, p, c, d, e, f, a)
    T = eltype(c)
    x = Vector{Tuple{Vector{T},Vector{T}}}()

    for j = 1:s-1
        push!(
            x,
            (vcat(zero(T), one(T), zeros(T, j - 1)), vcat(zeros(T, j), one(T))),
        )
    end

    # y0s
    push!(x, (vcat(zeros(T, s), one(T)), vcat(zero(T), c[1:s])))

    # first term y1s
    push!(
        x,
        (
            vcat(zero(T), d[1:s], one(T)),
            vcat(zero(T), zero(T), e[3:(s+1)], one(T)),
        ),
    )

    # y1s
    y1s = vcat(f[1:s+1], e[1], one(T))
    if (p == 0) || (p == s)
        z = y1s
        return graph_degopt(x, z)
    else #Apply PS-scheme as in (52)
        # Assumed to be an integer, according to paper, p = v*s
        v = convert(Int, p / s)

        # Evaluate y1s*x^s
        push!(x, (y1s, vcat(zeros(T, s), one(T), zeros(T, 2))))

        # Evaluate rest of the multiplications with PS-scheme
        for i in reverse(0:v-2)
            lower_idx = (i + 1) * s + 1
            upper_idx = lower_idx + s - 1
            idx = lower_idx:upper_idx
            c = view(a, idx)
            # Coeffs for P_{v-i} + A^s*P_{v+1-i}
            xl = vcat(c, zero(T), zeros(T, v - i), one(T))
            xr = vcat(zeros(T, s), one(T), zeros(T, v + 1 - i)) # Coeffs for A^s
            push!(x, (xl, xr)) # A^s*( P_{v-i} + A^s*P_{v+1-i} )
        end

        z = vcat(a[1:s], zeros(T, 2 + v), one(T))
        return graph_degopt(x, z)
    end
end

# Internal use only
# Transforms formula (34)-(35) + (69) to degopt form
# Input is analogous to graph_sastre_z1ps_degopt() but with tuples of vectors
# # Nof mult: s+4
function graph_sastre_h2m_degopt(s, c, d, e, f)
    T = eltype(eltype(c))
    x = Vector{Tuple{Vector{T},Vector{T}}}()

    for j = 1:s-1
        push!(
            x,
            (vcat(zero(T), one(T), zeros(T, j - 1)), vcat(zeros(T, j), one(T))),
        )
    end

    # y0s
    push!(x, (vcat(zeros(T, s), one(T)), vcat(zero(T), c[1][1:s])))
    # first term y1s
    push!(
        x,
        (
            vcat(zero(T), d[1][1:s], one(T)),
            vcat(zero(T), zero(T), e[1][3:(s+1)], one(T)),
        ),
    )
    # y1s (adjusted for two more multiplications)
    y1s = vcat(f[1][1:s+1], e[1][1], one(T), zeros(T, 2))

    # y0s_p
    push!(
        x,
        (
            vcat(zeros(T, s), one(T), zeros(T, 2)),
            vcat(zero(T), c[2][1:s], zeros(T, 2)),
        ),
    )
    # first term y1s_p
    push!(
        x,
        (
            vcat(zero(T), d[2][1:s], zeros(T, 2), one(T)),
            vcat(zero(T), zero(T), e[2][3:(s+1)], zeros(T, 2), one(T)),
        ),
    )
    # y1s_p (skipping two multiplications in the middle)
    y1s_p = vcat(f[2][1:s+1], zeros(T, 2), e[2][1], one(T))

    push!(x, (vcat(y1s), y1s_p))

    z = vcat(one(T), zeros(s + 4), one(T))
    return graph_degopt(x, z)
end
