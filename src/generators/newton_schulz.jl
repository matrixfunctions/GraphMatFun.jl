export graph_newton_schulz, graph_newton_schulz_degopt

"""
     (graph,crefs)=graph_newton_schulz(
        k,
        T = ComplexF64;
        input = :A,
        Z = :Z,
        Q = :Q,
        V = :V,
    )

Does `k` iterations of the Newton–Schulz iteration for approximating the inverse
of a matrix A (name given by `input`), i.e., the recursion

    V_k+1 = V_k*(2*I - A*V_k),

with `V_0=A`.
The recursion is implemented using the graph operations

    Z_i=A*V_i
    Q_i=2*I-Z_i
    V_{i+1}=V_i*Q_i,

and the kwargs `Z`, `Q`, and `V` determines the naming of the corresponding
intermediate steps.

**References**

1. G. Schulz. "Iterative Berechnung der reziproken Matrix". Zeitschrift für
   Angewandte Mathematik und Mechanik, 13:57–59, 1933.
   DOI: [10.1002/zamm.19330130111](https://doi.org/10.1002/zamm.19330130111)

2. N. J. Higham. "Functions of Matrices". SIAM, Philadelphia, PA, 2008.
   DOI: [10.1137/1.9780898717778](https://doi.org/10.1137/1.9780898717778)
"""
function graph_newton_schulz(
    k,
    T = ComplexF64;
    input = :A,
    Z = :Z,
    Q = :Q,
    V = :V,
)
    graph = Compgraph(T)
    cref = Vector{Tuple{Symbol,Int}}(undef, 0)

    for i = 1:k
        if (i == 1)
            Vi = input
        else
            Vi = Symbol("$(V)$i")
        end
        Zi = Symbol("$(Z)$i")
        add_mult!(graph, Zi, input, Vi)
        Qi = Symbol("$(Q)$i")
        add_lincomb!(graph, Qi, 2, :I, -1.0, Zi)
        Vii = Symbol("$(V)$(i+1)")
        add_mult!(graph, Vii, Vi, Qi)
        add_output!(graph, Vii)
    end
    Vii = Symbol("$(V)$(k+1)")

    return (graph, cref)
end

"""
     (graph,crefs)=graph_newton_schulz_degopt(k, T=ComplexF64; input=:A)

Does `k` iterations of the Newton–Schulz iteration for approximating the inverse,
using the recursion

    Z_i=A*V_i
    V_{i+1}=V_i*(2*I-Z_i).

The function makes a call to [`graph_degopt`](@ref), resulting in more
degrees of freedom in `crefs`. See also [`graph_newton_schulz`](@ref).
"""
function graph_newton_schulz_degopt(k, T = ComplexF64; input = :A)
    x = Vector{Tuple{Vector{T},Vector{T}}}(undef, 2 * k)
    # Z_i=A*V_i --- Odd numbers
    for i = 1:2:(2*k-1)
        x[i] =
            (vcat(zero(T), one(T), zeros(T, i - 1)), vcat(zeros(T, i), one(T)))
    end
    # V_{i+1}=V_i*(2*I-Z_i) --- Even numbers
    for i = 2:2:(2*k)
        x[i] = (
            vcat(zeros(T, i - 1), one(T), zero(T)),
            vcat(2, zeros(T, i - 1), -one(T)),
        )
    end

    return graph_degopt(x, vcat(zeros(T, 2 * k + 1), one(T)), input = input)
end
