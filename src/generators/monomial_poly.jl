export gen_monomial, gen_monomial_recursive

"""
     (graph,crefs)=gen_monomial(a; input=:A, polyname=:P)

Generates the graph for the polynomial using the monomial basis coefficients. More precisely,
it corresponds to the evaluation of the polynomial

    p(s) = a[1] + a[2]*s + ... + a[n]*s^(n-1),

where s^k is naively evaluated as s^k=s*s^(k-1) for k=2,3,...,n-1.
The kwarg `polyname` specifies the name of intermediate variables.
    """
function gen_monomial(a; input=:A, polyname=:P)

    # Initial setup
    n = length(a);
    cref=Vector{Tuple{Symbol,Int}}(undef,n)
    graph = Compgraph(eltype(a));
    d = n-1; # Polynomial degree

    # Degenerate case d = 0
    if d == 0
        outkey = Symbol("$(polyname)2")
        add_lincomb!(graph, outkey, a[1], :I, 0.0, :I)
        cref[1] = (outkey,1);
        add_output!(graph,outkey);
        return (graph,cref)
    end

    outkey = Symbol("$(polyname)$n")
    nodelist = Vector{Symbol}(undef,n)
    nodelist[1] = :I
    nodelist[2] = Symbol("$input")

    # Create monomial basis
    key = input
    for i=2:d
        prevkey = key
        key = Symbol("$input$i")
        add_mult!(graph, key, input, prevkey)
        nodelist[i+1] = key
    end

    # Sum up the polynomial
    cref[:] .= add_sum!(graph, outkey, a, nodelist, polyname)

    add_output!(graph,outkey);
    return (graph,cref)

end



"""
     (graph,crefs)=gen_monomial_recursive(a; input=:A)

Generates the same polynomial as `gen_monomial`, in the monomial basis.
However, it does so by wrapping a call to `gen_general_poly_recursion`, resulting in more
degrees of freedom in `crefs`.

    """
function gen_monomial_recursive(a; input=:A)

    n = length(a)
    d = n-1
    T = eltype(a)

    if d == 0
        error("Does not implement degree-zero polynomial.")
    end

    x = Vector{Tuple{Vector{T},Vector{T}}}(undef,d-1)
    for i = 1:d-1
        x[i] = ( vcat(zeros(T,i),one(T)), vcat(zero(T),one(T),zeros(T,i-1)) )
    end

    return gen_general_poly_recursion(x, a, input=input)

end
