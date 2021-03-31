export gen_horner

"""
     (graph,crefs)=gen_horner(a; input=:A, B=:B, C=:C, scaling=1.0)

Generates the graph for the polynomial using Horner's method. More precisely,
it corresponds to the evaluation of the polynomial

    p(s) = a[1] + a[2]*(αs) + ... + a[n-1]*(αs)^(n-2) + a[n]*(αs)^(n-1),

as the recursion

    Cj=s*B{j+1}
    Bj=a[j]*I+α*Cj,

where α=`scaling`.

The kwargs `B` and `C` specifies the base-names of these intermediate variables.
    """
function gen_horner(a; input=:A, B=:B, C=:C, scaling=1.0)

    # Initial setup
    n = length(a);
    cref=Vector{Tuple{Symbol,Int}}(undef,n)
    graph = Compgraph(eltype(a))
    out = Symbol("$(B)1")
    d = n-1 # Polynomial degree

    if d==0
        add_lincomb!(graph, out, a[1], :I, 0.0, :I)
        add_output!(graph, out)
        cref[1] = (out,1)
        return (graph,cref)
    end

    for j=d:-1:1
        Bj=Symbol("$B$j")
        Bjp1=Symbol("$B$(j+1)")
        Cj=Symbol("$C$j")
        if (j==d)
            add_lincomb!(graph, Bj, a[j], :I, a[j+1]*scaling, input)
            cref[j+1] = (Bj,1)
            cref[j] = (Bj,2)
        else
            #   Cj=X*B{j+1}
            #   Bj=a[j]*I+α*Cj
            add_mult!(graph, Cj, input, Bjp1)
            add_lincomb!(graph, Bj, a[j], :I, scaling, Cj)
            cref[j] = (Bj,1)
        end
    end
    push!(graph.outputs,Symbol("$(B)1"))
    return (graph,cref)

end
