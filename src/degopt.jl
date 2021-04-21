# Degree optimal polynomials

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
    k=size(degopt.x[end][1]);
    push!(degopt.x,(zeros(T,k+1), zeros(T,k+1)))
    push!(degopt.y,zeros(T,k+2));
end
