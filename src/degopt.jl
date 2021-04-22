# Degree optimal polynomials
export Degopt, grow!, get_degopt_crefs, scale!, square!;
struct Degopt{T}
    x::Vector{Tuple{Vector{T},Vector{T}}}
    y::Vector{T}
end


function Degopt(x,y)
    T=eltype(y)
    return Degopt{T}(x,y);
end

"""
    grow!(degopt::Degopt)

Increases the degopt by one multiplication.

"""
function grow!(degopt::Degopt{T}) where T
    k=size(degopt.x[end][1],1);
    @show T
    @show degopt.x
    @show typeof((zeros(T,k+1), zeros(T,k+1)))
    @show degopt.y

    println("Push!");
    push!(degopt.x,(zeros(T,k+1), zeros(T,k+1)))
    println("Push2");
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
