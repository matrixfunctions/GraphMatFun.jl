# Degree optimal polynomials
export Degopt, grow!, get_degopt_crefs, scale!, square!;

struct Degopt{T}
    x::Vector{Tuple{Vector{T},Vector{T}}}
    y::Vector{T}
end


"""
    degopt=Degopt(x,y)
    degopt=Degopt(graph)

    struct Degopt{T}
        x::Vector{Tuple{Vector{T},Vector{T}}}
        y::Vector{T}
    end


Creates an object representing a degree optimal polynomial, i.e., the coefficients in an evaluation scheme which maximizes the degree for a fixed number of multiplications. The object represents the coefficients in


    A0=I
    A1=A
    A2=(x *A0+x *A1)(x *A0+x *A1)
    A4=(x *A0+x *A1+x*A2)(x *A0+x *A1+x *A2)
    A8=(x *A0+x *A1+x*A2)(x *A0+x *A1+x *A2)
     ..

and

    Out=y*A0+y*A1+y*A4+y*A8+y*B4...

The `x`-values are given in the argument `x`, which is
a `Vector{Tuple{Vector{T},Vector{T}}}`, containing the elements
of each sum. The `z`-vector contains the elements
to form the output.

The coefficients in graphs generated by [`gen_degopt_poly`](@ref)
can be recovered by using `Degopt(graph)`.

"""
function Degopt(x,y)
    T=eltype(y)
    return Degopt{T}(x,y);
end

function Degopt(graph)
    # Number of multiplications
    k=count(values(graph.operations) .== :mult)
    # Check so it's not a complicated graph
    if (count(values(graph.operations) .== :ldiv)>0)
        error("The graph must originate from a Degopt form and cannot have divisions");
    end
    (x_cref,y_cref)=get_degopt_crefs(k)
    T=eltype(graph);
    xv=Vector{Tuple{Vector{T},Vector{T}}}(undef,k);
    y=Vector{T}(undef,k+2);

    for (i,x)=enumerate(x_cref)
        xv[i]=(get_coeffs(graph,x[1]),get_coeffs(graph,x[2]));
    end
    y=get_coeffs(graph,y_cref);
    return Degopt(xv,y)

end

"""
    grow!(degopt::Degopt)

Increases the degopt by one multiplication without modifying the function values.

"""
function grow!(degopt::Degopt{T}) where T
    k=size(degopt.x[end][1],1);
    push!(degopt.x,(zeros(T,k+1), zeros(T,k+1)))
    push!(degopt.y,zero(T));
end



"""
    (x,z)=get_degopt_crefs(k;compress_keys=false)

Retruns crefs related to `gen_degopt_poly`. Specifically
`x` is a `Vector{Tuple{Vector{Tuple{Symbol,Int}},Vector{Tuple{Symbol,Int}}}}`
such that `x[2][1]` corresponds to the coefficients of the left hand side of the
multiplication

    B2=(α_2_1 *I + α_2_2 *A)(β_2_1 *I + β_2_2 *A)

i.e., the crefs corresponding to `[α_2_1, α_2_2]`. See `gen_degopt_poly`. Hence,
`get_coeffs(graph,x[2][1])` returns the corresponding numerical values of the coefficients.

"""
function get_degopt_crefs(k;compress_keys=false)

    TT = Tuple{Symbol,Int}
    x=Vector{Tuple{Vector{TT},Vector{TT}}}(undef,k)
    for i = 1:k
        x[i] = ( Vector{TT}(undef,i+1), Vector{TT}(undef,i+1) )
    end

    for s=2:k+1

        for (b,base) = enumerate(["Ba$s", "Bb$s"])
            if s == 2
                x[s-1][b][1] = (Symbol(base), 1)
            end
            for i = 2:s-1
                if i == 2
                    x[s-1][b][1] = (Symbol("$(base)_2"), 1)
                end
                x[s-1][b][i] = (Symbol("$(base)_$(i)"), 2)
            end
            x[s-1][b][s] = (Symbol(base), 2)
        end
    end

    zlength = compress_keys ? 2 : k+2
    z = Vector{TT}(undef,zlength)
    z[1] = (:T2k2, 1)
    z[2] = (:T2k2, 2)
    for i = 3:zlength-1
        z[i] = ( Symbol("T2k$(i)"), 2 )
    end
    if !compress_keys
        z[k+2] = ( Symbol("T2k$(k+3)"), 2 )
    end

    return (x,z)
end


"""
    scale!(degopt::Degopt,α)

Effectively change a `Degopt` such that the input is scaled by `α`. If `p` is the original function, `p(α x)` will be the new function.

"""
function scale!(degopt::Degopt,α)
    for x=degopt.x
        x[1][2] *= α
        x[2][2] *= α
    end
    degopt.y[2] *= α;
end


"""
    square!(degopt::Degopt)

Effectively square a `Degopt` in the sense that the output is square. If `p` is the original function, `p(x)^2` will be the new function.

"""
function square!(degopt::Degopt)
    y0=degopt.y;
    push!(degopt.x,(deepcopy(y0),deepcopy(y0)));
    degopt.y[:]=zero(y0);
    push!(degopt.y,1);
end


function row1_normalize!(xx,a,b,c,d)
    q3=xx[3];
    xx[1] += q3*a*c;
    xx[2] += q3*(a*d+b*c);
    xx[3]=q3*b*d;
end
import LinearAlgebra.normalize!; # for overloading
"""
    normalize!(degopt,tp=:row1)

Normalizes the degopt coefficients, in the way specified by `tp`.
If the `rp==:row1` the degopt will be transformed
to an equivalent degopt with first row equal to `(0 1) (0 1)`.

"""
function normalize!(degopt::Degopt,tp=:row1)
    if (tp==:row1)

        a=degopt.x[1][1][1]
        b=degopt.x[1][1][2]
        c=degopt.x[1][2][1]
        d=degopt.x[1][2][2]

        for (i,x)=enumerate(degopt.x)
            if (i>1)
                for k=1:2
                    row1_normalize!(x[k],a,b,c,d);
                end
            end
        end
        row1_normalize!(degopt.y,a,b,c,d);


        degopt.x[1][1][1]=0
        degopt.x[1][1][2]=1;
        degopt.x[1][2][1]=0
        degopt.x[1][2][2]=1;

        return degopt;
    else
        error("Unknown normalization");
    end
end
