# Degree-optimal polynomials
export Degopt,
    grow!, get_degopt_crefs, scale!, normalize!, square!, get_degopt_coeffs;

struct Degopt{T}
    x::Vector{Tuple{Vector{T},Vector{T}}}
    y::Vector{T}
end

"""
    degopt=Degopt(x,y)
    degopt=Degopt(HA,HB,y)
    degopt=Degopt(graph)

    struct Degopt{T}
        x::Vector{Tuple{Vector{T},Vector{T}}}
        y::Vector{T}
    end

Creates an object representing a degree-optimal polynomial, i.e., the
coefficients in an evaluation scheme which maximizes the degree for a fixed
number of multiplications. The object represents the coefficients in

    A0=I
    A1=A
    A2=(x *A0+x *A1)(x *A0+x *A1)
    A4=(x *A0+x *A1+x*A2)(x *A0+x *A1+x *A2)
    A8=(x *A0+x *A1+x*A2)(x *A0+x *A1+x *A2)
     ..

and

    Out=y*A0+y*A1+y*A4+y*A8+y*B4...

The `x`-values are given in the argument `x`, which is a
`Vector{Tuple{Vector{T},Vector{T}}}`, containing the elements of each sum. The
`y`-vector contains the elements to form the output. In the constructor with
`HA`, `HB` and `y` the elements of `x` are stored as matrices.

The coefficients in graphs generated by [`graph_degopt`](@ref) can be recovered
by using `Degopt(graph)`.
"""
function Degopt(x, y)
    T = eltype(y)
    return Degopt{T}(x, y)
end

function Degopt(graph)
    # Number of multiplications
    k = count(values(graph.operations) .== :mult)
    # Check so it's not a complicated graph
    if (count(values(graph.operations) .== :ldiv) > 0)
        error(
            "The graph must originate from a Degopt form " *
            "and cannot have divisions",
        )
    end
    (x_cref, y_cref) = get_degopt_crefs(k)
    T = eltype(graph)
    xv = Vector{Tuple{Vector{T},Vector{T}}}(undef, k)
    y = Vector{T}(undef, k + 2)

    for (i, x) in enumerate(x_cref)
        xv[i] = (get_coeffs(graph, x[1]), get_coeffs(graph, x[2]))
    end
    y = get_coeffs(graph, y_cref)
    return Degopt(xv, y)
end

function Degopt(HA::AbstractMatrix, HB::AbstractMatrix, y::AbstractVector)
    m = size(HA, 1)
    T = eltype(HA)
    x = Vector{Tuple{Vector{T},Vector{T}}}(undef, m)
    for k = 1:m
        x[k] = (HA[k, 1:k+1], HB[k, 1:k+1])
    end
    return Degopt(x, y)
end

"""
    grow!(degopt::Degopt)

Increases the [`Degopt`](@ref) by one multiplication without modifying the
function values.
"""
function grow!(degopt::Degopt{T}) where {T}
    k = size(degopt.x[end][1], 1)
    push!(degopt.x, (zeros(T, k + 1), zeros(T, k + 1)))
    push!(degopt.y, zero(T))
    return nothing
end

"""
    (x,z)=get_degopt_crefs(k)
    (x,z)=get_degopt_crefs(graph)

Returns linear combination references (crefs) related
to [`graph_degopt`](@ref). Specifically
`x` is a `Vector{Tuple{Vector{Tuple{Symbol,Int}},Vector{Tuple{Symbol,Int}}}}`
such that `x[2][1]` corresponds to the coefficients of the left hand side of the
multiplication

    B2=(α_2_1 *I + α_2_2 *A)(β_2_1 *I + β_2_2 *A)

i.e., the crefs corresponding to `[α_2_1, α_2_2]`. See [`graph_degopt`](@ref).
Hence, `get_coeffs(graph,x[2][1])` returns the corresponding numerical values of
the coefficients.
"""
function get_degopt_crefs(k)
    TT = Tuple{Symbol,Int}
    x = Vector{Tuple{Vector{TT},Vector{TT}}}(undef, k)
    for i = 1:k
        x[i] = (Vector{TT}(undef, i + 1), Vector{TT}(undef, i + 1))
    end

    for s = 2:k+1
        for (b, base) in enumerate(["Ba$s", "Bb$s"])
            for i=1:s-1
                x[s-1][b][i]=(Symbol(base),i)
            end
        end
    end

    z = map(i-> Symbol(:y,i), 1:k+2)
    return (x, z)
end

function get_degopt_crefs(graph::Compgraph)
    return get_degopt_crefs(count(values(graph.operations) .== :mult))
end

"""
    scale!(degopt::Degopt,α)

Effectively change a [`Degopt`](@ref) such that the input is scaled by `α`. If
`p` is the original function, `p(α x)` will be the new function.
"""
function scale!(degopt::Degopt, α)
    for x in degopt.x
        x[1][2] *= α
        x[2][2] *= α
    end
    degopt.y[2] *= α
    return nothing
end

"""
    square!(degopt::Degopt)

Effectively square a [`Degopt`](@ref) in the sense that the output is square. If
`p` is the original function, `p(x)^2` will be the new function.
"""
function square!(degopt::Degopt)
    y0 = degopt.y
    push!(degopt.x, (deepcopy(y0), deepcopy(y0)))
    degopt.y[:] = zero(y0)
    push!(degopt.y, 1)
    return nothing
end

function row1_normalize!(xx, a, b, c, d)
    q3 = xx[3]
    xx[1] += q3 * a * c
    xx[2] += q3 * (a * d + b * c)
    xx[3] = q3 * b * d
    return nothing
end

import LinearAlgebra.normalize!; # for overloading
"""
    normalize!(degopt::Degopt,tp=:row1) -> degopt

Normalizes the [`Degopt`](@ref) coefficients, in the way specified by `tp`. If
the `tp==:row1` the `degopt` will be transformed to an equivalent
[`Degopt`](@ref) with first row equal to `(0 1) (0 1)`.  If `tp==:col1` the first column in the `Ha` and `Hb` matrices will be transformed to zero.
"""
function normalize!(degopt::Degopt, tp = :row1)
    if (tp == :row1)
        a = degopt.x[1][1][1]
        b = degopt.x[1][1][2]
        c = degopt.x[1][2][1]
        d = degopt.x[1][2][2]

        for (i, x) in enumerate(degopt.x)
            if (i > 1)
                for k = 1:2
                    row1_normalize!(x[k], a, b, c, d)
                end
            end
        end
        row1_normalize!(degopt.y, a, b, c, d)

        degopt.x[1][1][1] = 0
        degopt.x[1][1][2] = 1
        degopt.x[1][2][1] = 0
        degopt.x[1][2][2] = 1

        return degopt
    elseif (tp == :col1)
        (Ha,Hb,y)=get_degopt_coeffs(degopt);
        # For loop through all rows
        m=size(Ha,1);
        for j=1:m
            # For each row make the expansion
            # (ca*I+fa)(cb*I+fb)=ca*cb*I+cb*fa+ca*fb+ fa*fb
            # Let the new row be just fa*fb and compensate for
            #   ca*cb*I+cb*fa+ca*fb
            # in all subsequent rows in the Ha, Hb and y

            ca=Ha[j,1];
            cb=Hb[j,1];
            c=ca*cb;
            fa=[0;Ha[j,2:j+1]];
            fb=[0;Hb[j,2:j+1]];

            Ha[j,1]=0;
            Hb[j,1]=0;

            # The remainder/compansator to be inserted into later rows
            compensator=(cb*fa+ca*fb);
            compensator[1]+=c;
            for k=j+1:m
                Ha[k,1:j+1] += compensator*Ha[k,j+2];
                Hb[k,1:j+1] += compensator*Hb[k,j+2];
            end
            y[1:j+1] += compensator*y[j+2];
        end
        degopt2=Degopt(Ha,Hb,y);
        degopt.x[:]=degopt2.x[:];
        degopt.y[:]=degopt2.y[:];
        return degopt
    else
        error("Unknown normalization")
    end
end

"""
    get_degopt_coeffs(degopt) -> (HA,HB,y)
    get_degopt_coeffs(graph) -> (HA,HB,y)

Returns the coefficients of the degree-optimal polynomial using two matrices and
one vectors as described in [`Degopt`](@ref).
"""
function get_degopt_coeffs(degopt::Degopt)
    p = size(degopt.x, 1)
    T = eltype(degopt.y)
    HA = zeros(T, p, p + 1)
    HB = zeros(T, p, p + 1)
    for i = 1:p
        HA[i, 1:(i+1)] = degopt.x[i][1]
        HB[i, 1:(i+1)] = degopt.x[i][2]
    end
    y = degopt.y
    return (HA, HB, y)
end
get_degopt_coeffs(graph::Compgraph) = get_degopt_coeffs(Degopt(graph));
