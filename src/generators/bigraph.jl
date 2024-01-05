export graph_bigraph

# Function that helps us find a basename which is not already
# used in the graph
function unused_base_name(graph,inspired_by)

    ismodified = true
    while ismodified
        ismodified=false;
        if (any(map(n-> startswith(string(n),string(inspired_by)) && n != inspired_by, collect(keys(graph.parents)))))
            # problematic case
            inspired_by=Symbol(string(inspired_by)*"_bigraph");
            ismodified=true;
        end
    end
    return inspired_by;
end
"""
     (graph,cref)=graph_bigraph(graph)

Turns the graph into a graph where very node has exactly two parents.

"""
function graph_bigraph(graph0)
    graph=deepcopy(graph0);
    lincomb_keys=keys(graph.coeffs)
    for (i,node)=enumerate(lincomb_keys)
        coeffs=graph.coeffs[node];
        parents=graph.parents[node];

        if length(coeffs)>2
            # Overwrite it!
            add_sum!(graph,node, coeffs,parents,
                     unused_base_name(graph,node))
        elseif (length(coeffs)==1) # Expand by a trivial multiplication
            graph.parents[node]=[:I;parents[1]];
            graph.coeffs[node]=[0;coeffs[1]];
        end
    end
    cref=get_all_cref(graph);
    return (graph,cref)

end
