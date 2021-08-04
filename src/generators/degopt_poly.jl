export graph_degopt, get_topo_order_degopt

"""
    (graph,crefs)=graph_degopt(k;T=ComplexF64,input=:A)
    (graph,crefs)=graph_degopt(x,z;input=:A)
    (graph,crefs)=graph_degopt(d::Degopt;input=:A)

Corresponds to the (for a fixed numer of multiplications) degree-optimal
polynomial

    B1=A
    B2=(x *I+x *A)(x *I+x *A)
    B3=(x *I+x *A+x*B2)(x *I+x *A+x *B2)
     ..

and

    Out=z*I+z*A1+z*B2+z*B3+z*B4...

The `x`-values are given in the argument `x`, which is a
`Vector{Tuple{Vector,Vector}}`, containing the elements of each sum. The
`z`-vector contains the elements to form the output, and `input` determines the
name of the matrix A above. If the parameter `k` is supplied instead of the
coefficients, all coeffs will be set to one. The general recursion
[Equation (9), 1].

Reference:

[1] P. Bader, S. Blanes, and F. Casas, "Computing the matrix exponential
    with an optimized Taylor polynomial approximation", Mathematics, 7(12),
    2019. DOI: [10.3390/math7121174](https://doi.org/10.3390/math7121174)
"""
function graph_degopt(x, z; input = :A)
    T = promote_type(eltype(eltype(eltype(x))), eltype(z))
    (graph, crefs) = graph_degopt_B(x, T, input = input)

    k = size(x, 1)

    # Add the polynomial in the Bk coeffs

    z_nodes = [:I; input]
    for s = 2:k+1
        push!(z_nodes, Symbol("B$s"))
    end
    key = Symbol("T2k$(k+3)")
    crefs_new = add_sum!(graph, key, z, z_nodes, Symbol("T2k"))
    append!(crefs, crefs_new)

    # Set the output
    empty!(graph.outputs)
    push!(graph.outputs, key)

    return (graph, crefs)
end
function graph_degopt(k; T = ComplexF64, input = :A)
    x = Vector{Tuple{Vector{T},Vector{T}}}(undef, k)
    for i = 1:k
        x[i] = (ones(T, i + 1), ones(T, i + 1))
    end

    z = ones(T, k + 2)

    return graph_degopt(x, z; input = :A)
end
function graph_degopt(degopt::Degopt; input = :A)
    return graph_degopt(degopt.x, degopt.y, input = input)
end

# Normally x::Vector{Tuple{Vector{Number},Vector{Number}}}
# Containing the coefficients in the recursion
function graph_degopt_B(x, T; input = :A)
    k = size(x, 1)

    graph = Compgraph(T)
    useful_syms = [:I, input]

    cref = []
    for s = 2:k+1
        base0 = "B$s"
        base_a = "Ba$s"
        base_b = "Bb$s"

        # First poly
        c_a = x[s-1][1]
        crefs_a = add_sum!(
            graph,
            Symbol(base_a),
            c_a,
            useful_syms,
            Symbol("$(base_a)_"),
        )

        append!(cref, crefs_a)

        # Second poly
        c_b = x[s-1][2]
        crefs_b = add_sum!(
            graph,
            Symbol(base_b),
            c_b,
            useful_syms,
            Symbol("$(base_b)_"),
        )
        append!(cref, crefs_b)

        # Multiply them together
        add_mult!(graph, Symbol(base0), Symbol(base_a), Symbol(base_b))
        push!(useful_syms, Symbol(base0))
        empty!(graph.outputs)
        push!(graph.outputs, Symbol(base0))
    end

    return (graph, cref)
end

"""
    order=get_topo_order_degopt(k)

A special implementation of [`get_topo_order`](@ref) for degree-optimal polynomials
generated with [`graph_degopt`](@ref). The natural order of computation is to compute
row by row.

See also [`get_degopt_crefs`](@ref).
"""
function get_topo_order_degopt(k)
    (x, z) = get_degopt_crefs(k)
    computation_order = Vector{Symbol}(undef, 0)
    for i = 1:k
        for n = 1:2
            for j = 2:(i+1)
                push!(computation_order, x[i][n][j][1])
            end
        end
        push!(computation_order, Symbol("B$(i+1)"))
    end
    for i = 2:(2+k)
        push!(computation_order, z[i][1])
    end
    return computation_order
end
