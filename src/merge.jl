export merge_graphs;

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

    return Compgraph(operations, parents, coeffs, outputs)
end
