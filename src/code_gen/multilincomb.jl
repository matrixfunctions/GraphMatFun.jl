# Represents a graph which with different lincomb functionality:
# It can add more than two matrices. Only used in code generation.
struct MultiLincombCompgraph{T}
    operations::Dict{Symbol,Symbol}
    parents::Dict{Symbol,NTuple{<:Any,Symbol}}
    coeffs::Dict{Symbol,NTuple{<:Any,T}}
    outputs::Vector{Symbol}
end

function MultiLincombCompgraph(g::Compgraph)
    T = eltype(g)
    extgraph = MultiLincombCompgraph(
        Dict{Symbol,Symbol}(),
        Dict{Symbol,NTuple{<:Any,Symbol}}(),
        Dict{Symbol,NTuple{<:Any,T}}(),
        Vector{Symbol}(),
    )
    Z = extract_sums(g)
    for k in keys(g.operations)
        if (g.operations[k] == :mult)
            add_mult!(extgraph, k, g.parents[k][1], g.parents[k][2])
        end
        if (g.operations[k] == :ldiv)
            add_ldiv!(extgraph, k, g.parents[k][1], g.parents[k][2])
        end
    end
    for s in Z
        coeff_list = s[1]
        symbol_list = s[2]
        key = s[3][end]
        p = size(coeff_list, 1)
        extgraph.operations[key] = :lincomb
        extgraph.coeffs[key] = NTuple{p,T}(coeff_list)
        extgraph.parents[key] = NTuple{p,Symbol}(symbol_list)
    end
    for k in g.outputs
        push!(extgraph.outputs, k)
    end
    return extgraph
end


"""
    sums=extract_sums(graph::Compgraph{T})

`sums::Vector{Tuple{Vector{T},Vector{Symbol},Vector{Symbol}}}`

Returns a representation of sums in `graph` which may potentially be merged, in a dot-fusion.
The vector `sums` contains a tuple for each of these sums. The three entries of
the tuple are:

  - a vector of `T` values that represent the coefficients of the summands;
  - a vector of `Symbol`s that correspond to the summands; and
  - a vector of intermediate `Symbol`s (i.e, nodes) that can be merged.

The first two vectors have the same number of entries, one for each element that
can be merged in the sum.
"""
function extract_sums(graph)
    coeff, nodes, merged, sums =
        find_mergeable_sums(graph, graph.outputs[1], [])
    return sums
end

# Return true if `node` has multiple parents. Such nodes cannot be freed.
function has_multiple_parents(graph, node)
    return sum(map(x -> any(x .== node), values(graph.parents))) > 1
end

function find_mergeable_sums(graph, node, processed, curr_coeff = 1)
    # Extract sums in the subgraph of `graph` with root `node`.
    #
    # The function accepts four parameters:
    #     * `graph`: current graph.
    #     * `node`: current node.
    #     * `processed`: set of nodes the algorithm has already processed..
    #     * `curr_coeff`: coefficient `node` if parent is lincomb, 1 otherwise.
    #
    # The functions returns four vectors:
    #     * `pcoeffs`: coefficients of sum currently being constructed.
    #     * `pnodes`: corresponding nodes of sum currently being constructed.
    #     * `pmerged`: nodes merged in the current sum (these will disappear).
    #     * `sums`: completely extracted sums.
    #
    # The algorithm starts from the first output node, which is seen as the root
    # of a spanning tree with edges defined in `graph.parents`. The graph may
    # have cycles, but the vector `processed` ensures that each node is
    # processed only once, the first time it is visited.
    #
    # If `node` is an input node, that is, a node without parents, then the
    # algorithm returns four empty vectors, as 1) the subtree rooted at `node`,
    # being empty, does not have extracted sums, and 2) no sum is being
    # constructed.
    #
    # If the node is not a leaf, the function is called recursively on the
    # two parents, and three cases are possible:
    #
    # 1) If `node` is not a `:lincomb`, then the function returns the union of
    # the sums extracted in the two subgraphs rooted at the parents. If either
    # parent is a `:lincomb` the sum that parent was constructing is added to
    # the vector of extracted sums. The three other output vectors are empty.
    #
    # 2) If `node` is a `:lincomb`, is parent to only one node, and is not an
    # output node, then the function merges the two sums being constructed by
    # the parents, if any, adds `node` to it, and returns the data accordingly.
    # The union of the vectors of extracted sums is also returned.
    #
    # 3) Otherwise, `node` is added to the union of the (possibly empty) sums
    # being constructed by the parents. In particular, the algorithm will add
    # `node` to the list of nodes to be merged, will updated the coefficients of
    # the constructed sum accordingly, and will add the current sum to the
    # vector of extracted sums, which will also include the union of the sum
    # extracted in the subgraph rooted at the parents.

    if !(node in keys(graph.operations)) || (node in processed)
        # Nothing to do for leaf nodes and nodes already processed.
        return Float64[], Symbol[], Symbol[], []
    else
        # Call function recursively on both parents.
        push!(processed, node)
        (parent1, parent2) = graph.parents[node]
        curr_lincomb = graph.operations[node] == :lincomb
        (coeff1, coeff2) = curr_lincomb ? graph.coeffs[node] : (1, 1)
        pcoeffs1, pnodes1, pmerged1, sums1 =
            find_mergeable_sums(graph, parent1, processed, coeff1)
        pcoeffs2, pnodes2, pmerged2, sums2 =
            find_mergeable_sums(graph, parent2, processed, coeff2)
        if curr_lincomb
            # Grow the sum by adjoining terms coming from parents, if any.
            has_multiple_parents(graph, node) && (curr_coeff = 1)
            new_coeffs = vcat(curr_coeff * pcoeffs1, curr_coeff * pcoeffs2)
            new_nodes = vcat(pnodes1, pnodes2)
            new_merged = vcat(pmerged1, pmerged2)
            # Lincomb parents are added to the mergeable nodes.
            # Non-lincomb parents are added to the sum and put, and their
            # coefficients are added to the vector of coefficients.
            if haskey(graph.operations, parent1) &&
               graph.operations[parent1] == :lincomb &&
               !has_multiple_parents(graph, parent1)
                new_merged = vcat(new_merged, parent1)
            else
                new_coeffs = vcat(new_coeffs, curr_coeff * coeff1)
                new_nodes = vcat(new_nodes, parent1)
            end
            if haskey(graph.operations, parent2) &&
               graph.operations[parent2] == :lincomb &&
               !has_multiple_parents(graph, parent2)
                new_merged = vcat(new_merged, parent2)
            else
                new_coeffs = vcat(new_coeffs, curr_coeff * coeff2)
                new_nodes = vcat(new_nodes, parent2)
            end
            if node in graph.outputs || has_multiple_parents(graph, node)
                return Float64[],
                Symbol[],
                Symbol[],
                vcat(
                    sums1,
                    sums2,
                    (new_coeffs, new_nodes, vcat(new_merged, node)),
                )
            else
                return new_coeffs, new_nodes, new_merged, vcat(sums1, sums2)
            end
        else
            # Current node is not a lincomb node.
            sums = vcat(sums1, sums2)
            # Add sum of lincomb parents, if any, to sums.
            sums =
                !isempty(pcoeffs1) ?
                vcat(sums, (pcoeffs1, pnodes1, vcat(pmerged1, parent1))) : sums
            sums =
                !isempty(pcoeffs2) ?
                vcat(sums, (pcoeffs2, pnodes2, vcat(pmerged2, parent2))) : sums
            return Float64[], Symbol[], Symbol[], sums
        end
    end
end
