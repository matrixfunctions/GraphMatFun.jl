export graph_degopt

"""
    (graph,crefs)=graph_degopt(k;T=ComplexF64,input=:A)
    (graph,crefs)=graph_degopt(x,z;input=:A)
    (graph,crefs)=graph_degopt(d::Degopt;input=:A)

Corresponds to the (for a fixed number of multiplications) degree-optimal
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
coefficients, all coeffs will be set to one. The general recursion is given in
(9) in the paper referenced below.

**Reference**

1. P. Bader, S. Blanes, and F. Casas, "Computing the matrix exponential
   with an optimized Taylor polynomial approximation", Mathematics, 7(12), 2019.
   DOI: [10.3390/math7121174](https://doi.org/10.3390/math7121174)
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
    key = Symbol("y")
    crefs_new = add_lincomb!(graph, key, z, z_nodes)
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
        crefs_a = add_lincomb!(
            graph,
            Symbol(base_a),
            c_a,
            useful_syms)

        append!(cref, crefs_a)

        # Second poly
        c_b = x[s-1][2]
        crefs_b = add_lincomb!(
            graph,
            Symbol(base_b),
            c_b,
            useful_syms)
        append!(cref, crefs_b)

        # Multiply them together
        add_mult!(graph, Symbol(base0), Symbol(base_a), Symbol(base_b))
        push!(useful_syms, Symbol(base0))
        empty!(graph.outputs)
        push!(graph.outputs, Symbol(base0))
    end

    return (graph, cref)
end
