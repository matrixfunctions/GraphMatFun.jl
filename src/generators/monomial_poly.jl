export gen_monomial

"""
     (graph,crefs)=gen_monomial(a; input=:A, polyname=:P)

Generates the graph for the polynomial using the monomial basis coefficieents. More precisely, 
it corresponds to evaluation of the polynomial

    p(s)=a[1]*I+a[2]*A+...a[n]*A^(n-1),

where A^k is naively evaluated as A^k=A*A^(k-1) for k=2,3,...,n-1. The kwarg `polyname` specifies the name of intermediate variables.
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
    nodelist[2] = :A

    # Create monomial basis
    A = input;
    key = input
    for i=2:d
        prevkey = key
        key = Symbol("$input$i")
        add_mult!(graph, key, A, prevkey)
        nodelist[i+1] = key
    end

    # Sum up the polynomial
    cref[:] .= add_sum!(graph, outkey, a, nodelist, polyname)

    add_output!(graph,outkey);
    return (graph,cref)

end
