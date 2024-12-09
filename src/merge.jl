export merge_graphs;
export split_lincomb!;

# Change a symbol name using prefix & postfix
function translate_keyname(key, prefix, postfix, skip_basic, input)
    newkey = Symbol(prefix * string(key) * postfix)
    if (skip_basic & (key == :I || key == input))
        newkey = key
    end
    return newkey
end

function change_keywords!(
    graph;
    prefix = "G",
    postfix = "",
    skip_basic = true,
    input = :A,
)
    translate =
        key -> translate_keyname(key, prefix, postfix, skip_basic, input)

    # Go through all dicts and update symbols according to translate
    for dict_org in (graph.coeffs, graph.operations, graph.parents)
        dict = deepcopy(dict_org)
        for (key, value) in dict
            delete!(dict_org, key)
            new_key = translate(key)
            dict_org[new_key] = value
            if (dict_org == graph.parents) # Also translate parents
                # Value is a Vector of symbols
                val=translate.(value)
                dict_org[new_key] = collect(val)
            end
        end
    end
    # Update outputs
    return graph.outputs[:] = translate.(graph.outputs)
end

"""
    graph=merge_graphs(
    graph1,
    graph2;
    prefix1 = "",
    prefix2 = "G2",
    skip_basic1 = true,
    skip_basic2 = true,
    cref1 = Vector(),
    cref2 = Vector(),
    input1 = :A,
    input2 = :A,

)

Takes all the nodes and edges in `graph1` and `graph2` and generates a new
graph. The node names are in `graph1` are changed by adding a prefix `prefix1`
and `graph2` correspondingly. The nodes `:I` and `:A` are unchanged if
`skip_basicX=true`. All coefficient references `crefX` are modified accordingly.
"""
function merge_graphs(
    graph1,
    graph2;
    prefix1 = "",
    prefix2 = "G2",
    skip_basic1 = true,
    skip_basic2 = true,
    cref1 = Vector(),
    cref2 = Vector(),
    input1 = :A,
    input2 = :A,
)
    T1 = eltype(graph1)
    T2 = eltype(graph2)

    T = promote_type(T1, T2)

    g1 = deepcopy(graph1)
    change_keywords!(
        g1,
        prefix = prefix1,
        skip_basic = skip_basic1,
        input = input1,
    )
    g2 = deepcopy(graph2)
    change_keywords!(
        g2,
        prefix = prefix2,
        skip_basic = skip_basic2,
        input = input2,
    )

    operations = merge(g1.operations, g2.operations)
    parents = merge(g1.parents, g2.parents)

    # Make the coeff types the same
    g1=Compgraph(T,g1);
    g2=Compgraph(T,g2);
    coeffs = merge(g1.coeffs, g2.coeffs);

    outputs = vcat(g1.outputs, g2.outputs)

    cref1_copy = copy(cref1)
    cref2_copy = copy(cref2)

    empty!(cref1)
    empty!(cref2)


    for c in cref1_copy
        push!(cref1, (Symbol(prefix1 * String(c[1])), c[2]))
    end
    for c in cref2_copy
        push!(cref2, (Symbol(prefix2 * String(c[1])), c[2]))
    end

    return Compgraph{T}(operations, parents, coeffs, outputs)
end


"""
      (g,cref_modified,new_crefs)=split_lincomb!(g,node,ind2;
                        newnode=node_new,
                        cref_list=[])

Takes the lincomb operation in node and splits it into two lincombs, by
creating a new node newnode. The new node consists of the linear combination
of the coefficients moving (node, ind2[1]),... (node,ind2[end]) and the
old lincomb object has the old coefficients and an additional term pointing
to newnode. The cref_list is updated to the new coefficient pointers. The
new_crefs list contains all the new coefficient pointers.


In this way, the graph is unchanged but one of the linear combinations is
split up into two.

"""
function split_lincomb!(g,node,ind2;
                        newnode=Symbol("$(node)_new"),
                        cref_list=[])

    @show newnode
    # get ind1 = complement of ind2
    nof_lincombs=size(g.coeffs[node],1);
    @show nof_lincombs
    ind1=map(i -> !(i in ind2), 1:nof_lincombs)

    org_parents=g.parents[node];
    org_coeffs=g.coeffs[node];

    @show org_coeffs

    # Move ind2 to a new lincomb
    add_lincomb!(g,newnode,org_coeffs[ind2],org_parents[ind2])

    @show newnode
    # Store ind1 lincomb info
    new_parents1=[org_parents[ind1];newnode]
    new_coeffs1=[org_coeffs[ind1];1]

    # Update the node lincomb data to point ind2 + newnode
    empty!(g.coeffs[node]);
    push!(g.coeffs[node], new_coeffs1...)
    empty!(g.parents[node]);
    push!(g.parents[node], new_parents1...)


    # Update the cref_list

    new_crefs=[];
    replace_list=Dict();
    for (j,i)=enumerate(ind2);
        replace_list[(node,i)]=(newnode,j)
        push!(new_crefs,(newnode,j));
    end


    for (cref_old,cref_new) in replace_list
        @show cref_old
        @show cref_new
        ii=findall( [cref_old] .== cref_list )
        map(j -> cref_list[j]=cref_new, ii);
    end

    return (g,cref_list,new_crefs)
end
