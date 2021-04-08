export gen_newton_schulz

"""
     (graph,crefs)=gen_newton_schulz(k, T=ComplexF64; input=:A, B=:B, C=:C, V=:V)

Does `k` iterations of the Newton–Schulz iteration for approximating the inverse
of a matrix A (name given by `input`), i.e., the recursion

    V_k+1 = V_k*(2*I - A*V_k),

with V_0=A.
The recursion is implemented using the graph-operations

    Z_i=A*V_i
    Q_i=2*I-Z_i
    V_{i+1}=V_i*Q_i,

and the kwargs  `Z`, `Q`, and `V` determines the naming of the corresponding
intermediate steps.

References:

Günther Schulz. Iterative Berechnung der reziproken Matrix. Z. Angew. Math. Mech.,
13:57–59, 1933.

N. J. Higham. Functions of Matrices. SIAM publications, Philadelphia, PA, 2008.
    """
function gen_newton_schulz(k, T=ComplexF64; input=:A, Z=:Z, Q=:Q, V=:V)

    graph = Compgraph(T);
    cref=Vector{Tuple{Symbol,Int}}(undef,0)

    for i=1:k

        if (i==1)
            Vi = input
        else
            Vi = Symbol("$(V)$i")
        end
        Zi = Symbol("$(Z)$i")
        add_mult!(graph,Zi, input, Vi);
        Qi = Symbol("$(Q)$i")
        add_lincomb!(graph, Qi, 2, :I, -1.0, Zi)
        Vii = Symbol("$(V)$(i+1)")
        add_mult!(graph, Vii, Vi, Qi);
        add_output!(graph,Vii);
    end
    Vii = Symbol("$(V)$(k+1)")

    return (graph,cref)
end
