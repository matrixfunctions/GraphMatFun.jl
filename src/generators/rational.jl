export gen_rational
"""
     (graph, cref) = gen_rational(den_coeffs, num_coeffs, poly_gen=gen_ps)
     (graph, cref) = gen_rational(den_graph, den_graph; den_cref=Vector{Tuple{Symbol,Int}}(), num_cref=Vector{Tuple{Symbol,Int}}())

Generates the graph for the rational approximation

    r(A)=q(A)^{-1}p(A)

where p(A) and q(A) are polynomials defined by the coeficients `den_coeffs` and `num_coeffs`,
and generated by the function `poly_gen`, which is called as (graph,cref)=poly_gen(coeffs),
see, e.g., `gen_monomial` and `gen_ps`.

The alternative call-signature involves the graphs for p and q directly, as `den_graph`
and `num_graph`. The corresponding `den_cref` and `num_cref` can also be passed
to be modified accodringly, otherwise the return value `cref` is empty for this call.
    """
function gen_rational(den_coeffs, num_coeffs, poly_gen)

    (den_graph, den_cref) = poly_gen(den_coeffs)
    (num_graph, num_cref) = poly_gen(num_coeffs)

    return gen_rational(den_graph, num_graph; den_cref=den_cref, num_cref=num_cref)
end

function gen_rational(den_graph, num_graph; den_cref=Vector{Tuple{Symbol,Int}}(), num_cref=Vector{Tuple{Symbol,Int}}())

    graph = merge_graphs(den_graph, num_graph, prefix1="de", prefix2="nu", skip_basic1=true, skip_basic2=true, cref1=den_cref, cref2=num_cref)
    add_ldiv!(graph, :pade, graph.outputs[2], graph.outputs[1])
    cref =vcat(den_cref, num_cref);
    empty!(graph.outputs)
    push!(graph.outputs, :pade)

    return (graph, cref)
end