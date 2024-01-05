# Functions to automatically remove
# nodes from graphs

export compress_graph_dangling!
export compress_graph_zero_coeff!
export compress_graph_output_cleaning!
export has_identity_lincomb
export has_trivial_nodes
export compress_graph_trivial!
export compress_graph_redundant!
export compress_graph_passthrough!
export compress_graph!

function conditional_println(out_string, verbose)
    return verbose && println(out_string)
end

function delete_crefs!(cref, n)
    II = [] # II will be list of nodes not to remove
    for (i, c) in enumerate(cref)
        if (c[1] != n)
            push!(II, i)
        end
    end
    cref_org = copy(cref)
    empty!(cref)
    return append!(cref, cref_org[II])
end

"""
    compress_graph_dangling!(graph,cref=[];verbose=false)

Removes dangling nodes in the graph. If a node is not used anywhere and it is
not in the output, it can safely be removed.
"""
function compress_graph_dangling!(graph, cref = []; verbose = false)
    # Remove dangling nodes
    modified = true
    while modified
        modified = false
        for (key, parents) in graph.parents
            nof_children = 0
            for (pkey, pparents) in graph.parents
                nof_children += sum(pparents .== key)
            end

            if (nof_children == 0 && !(key in graph.outputs))
                conditional_println("Delete dangling: $key", verbose)
                del_node!(graph,key)
                delete_crefs!(cref, key)
                modified = true
                break
            end
        end
    end
end

"""
    compress_graph_zero_coeff!(graph,cref=[];droptol=0;verbose=false)

Searches for linear combinations with zero coeff and removes those. The cref
list deletes references to zero coeffs and updates all other crefs.
"""
function compress_graph_zero_coeff!(
    graph,
    cref = [];
    droptol = 0,
    verbose = false,
)
    modified = true
    while (modified)
        modified = false
        all_crefs=get_all_cref(graph);
        vals=get_coeffs(graph,  all_crefs)
        i=findfirst(abs.(vals) .<= droptol )
        if (!isnothing(i)) # We found a zero element
            node=all_crefs[i][1];
            idx=all_crefs[i][2];

            conditional_println(
                "Remove node: One coefficient in the computation of $node can be removed " *
                string(vals[i]) *
                " ≈ 0",
                verbose,
            )

            deleteat!(graph.parents[node],idx);
            deleteat!(graph.coeffs[node],idx);
            modified=true;


            # Update the cref list

            # Remove the occurance of this cref
            rm_cref=all_crefs[i];
            filter!(x->x≠rm_cref,cref)
            # Update all the pointers in cref list
            for j = findall(map(x-> x[1] == node, cref))
                if (cref[j][2] > idx)
                    cref[j] = (cref[j][1],cref[j][2]-1) # Reduce by one. Tuple is immutable
                end
            end
            #
        end
    end
end

# Keep in order to not lose functionality.
function compress_graph_zero_coeff_old!(
    graph,
    cref = [];
    droptol = 0,
    verbose = false,
)
    modified = true
    while (modified)
        modified = false
        for (key, parents) in graph.parents
            big_break = false
            if (graph.operations[key] == :lincomb && !(key in graph.outputs)) # don't remove outputs
                v = graph.coeffs[key]
                for s = 1:2
                    other_s = 3 - s
                    if (abs(v[s]) <= droptol)
                        conditional_println(
                            "Remove node: $key because coeff $s is " *
                            string(v[s]) *
                            " ≈ 0",
                            verbose,
                        )

                        # Redirect all lincombs we can find
                        for (pkey, op) in graph.operations
                            pp = [graph.parents[pkey]...]
                            if (
                                op == :lincomb &&
                                (pp[1] == key || pp[2] == key)
                            )
                                # It can be merged

                                conditional_println(
                                    "Delete connection: $key -> $pkey",
                                    verbose,
                                )
                                # z: parent connected to the node to be removed
                                z = 1
                                if (pp[2] == key)
                                    z = 2
                                end

                                cc = [graph.coeffs[pkey]...]
                                cc[z] = cc[z] * v[other_s]

                                pp[z] = graph.parents[key][other_s]
                                # Overwrite
                                add_lincomb!(
                                    graph,
                                    pkey,
                                    cc[1],
                                    pp[1],
                                    cc[2],
                                    pp[2],
                                )

                                modified = true
                                big_break = true
                            end
                        end
                        if (big_break)
                            break
                        end
                    end
                end
            end
            if (big_break)
                break
            end
        end
    end
end
"""
    compress_graph_output_cleaning!(graph,cref=[];verbose=false)

Checks if the output is computed from the linear combination that can be
compressed.
"""
function compress_graph_output_cleaning!(graph, cref = []; verbose = false)
    ismodified = false
    while ismodified
        ismodified = false
        for (i, n) in enumerate(graph.outputs)
            if (graph.operations[n] == :lincomb)
                for s = 1:2
                    other_s = 3 - s
                    if (
                        graph.coeffs[n][s] == 1.0 &&
                        graph.coeffs[n][other_s] == 0.0
                    )
                        # Redirect since it is just a copy of
                        # another node
                        graph.outputs[i] = graph.parents[n][s]
                        ismodified = true
                    end
                end
            end
        end
    end
end

### Remove trivial nodes:
# 1) :mult with identity parent;
# 2) :ldiv with identity as first (left) parent.

# Check if node is trivial because of an identity on the left (:mult or :ldiv).
function is_trivial_left(graph, node)
    op = graph.operations[node]
    if (op == :lincomb)
        return false;
    else
        return graph.parents[node][1] == :I &&
            (op == :mult || op == :ldiv)
    end
end

# Check if node is trivial because of an identity on the right (:mult).
function is_trivial_right(graph, node)

    op = graph.operations[node]
    if (op == :lincomb)
        return false;
    else
        return (graph.parents[node][2] == :I && op == :mult)
    end
end

# Return a parent that can replace a node.
# The empty array Any[] is returned if the node is not trivial.
function find_replacement_node(graph, node)
    if is_trivial_left(graph, node)
        return graph.parents[node][2]
    elseif is_trivial_right(graph, node)
        return graph.parents[node][1]
    else
        return []
    end
end

# Replace node with replacement node.
function replace_node!(graph, node, replacement, cref)
    # Make replacement an output node if current node was.
    if (node in graph.outputs)
        deleteat!(graph.outputs, graph.outputs .== node)
        add_output!(graph, replacement)
    end
    # Update children of node to point to replacement.
    for child in get_children(graph, node)
        graph.parents[child] =
            map(x -> (x == node ? replacement : x), graph.parents[child])
    end
    # Delete node from graph.
    del_node!(graph, node)
    return delete_crefs!(cref, node)
end

"""
    has_identity_lincomb(graph) -> Bool

Checks whether the graph has a node which is a linear combination of identity
matrices.
"""
function has_identity_lincomb(graph)
    has_identity_lincomb = false
    for (key, parents) in graph.parents
        if graph.operations[key] == :lincomb &&
           all(parents .== :I)
            has_identity_lincomb = true
            break
        end
    end
    return has_identity_lincomb
end

"""
    has_trivial_node(graph) -> Bool

Checks whether the graph has trivial nodes, that is, multiplications by the
identity or linear systems whose coefficient is the identity matrix.
"""
function has_trivial_nodes(graph)
    has_trivial_nodes = false
    for (key, parents) in graph.parents
        if is_trivial_left(graph, key) || is_trivial_right(graph, key)
            has_trivial_nodes = true
            break
        end
    end
    return has_trivial_nodes
end

"""
    compress_graph_trivial!(graph,cref=[];verbose=false)

Removes from graph the following operations
   I\\B -> B
   I*B -> B
   B*I -> B
"""
function compress_graph_trivial!(graph, cref = []; verbose = false)
    ismodified = true
    while ismodified
        ismodified = false
        for (key, parents) in graph.parents
            replacement = find_replacement_node(graph, key)
            if replacement != []
                conditional_println(
                    "Replace node: $key by $replacement",
                    verbose,
                )
                replace_node!(graph, key, replacement, cref)
                ismodified = true
            end
        end
    end
end

### Remove redundant nodes:
# 1) :mult with same coefficients;
# 2) :ldiv with same terms;
# 3) :lincomb with same terms and coefficients (not necessarily)

# Return list of nodes that can be merged with node.
function find_redundant_nodes(graph, node, compress_lincomb)
    # Find nodes with same operation and same parents.
    operation = graph.operations[node]
    if !compress_lincomb && operation == :lincomb
        return Array{Symbol}([])
    end
    parents = graph.parents[node]
    same_operation = keys(
        filter(
            kv -> kv.first != node && operation == kv.second,
            graph.operations,
        ),
    )
    if operation == :ldiv # Order of parents matter.
        redundant_candidates = keys(
            filter(
                kv -> kv.first in same_operation && kv.second == parents,
                graph.parents,
            ),
        )
    elseif operation == :mult # Order of parents doesn't matter.
        redundant_candidates = keys(
            filter(
                kv ->
                    kv.first in same_operation && (
                        kv.second == parents || kv.second == reverse(parents)
                    ),
                graph.parents,
            ),
        )
    elseif operation == :lincomb
        # For :lincomb nodes, check that coefficients are the same.
        coefficients = graph.coeffs[node]
        # Parents in the same order.
        same_parents_order = keys(
            filter(
                kv -> kv.first in same_operation && parents == kv.second,
                graph.parents,
            ),
        )
        redundant_candidates_order = keys(
            filter(
                kv ->
                    kv.first in same_parents_order &&
                        kv.second == coefficients,
                graph.coeffs,
            ),
        )
        # Parents in reverse order.
        same_parents_reverse = keys(
            filter(
                kv ->
                    kv.first in same_operation &&
                        parents == reverse(kv.second),
                graph.parents,
            ),
        )
        redundant_candidates_reverse = keys(
            filter(
                kv ->
                    kv.first in same_parents_reverse &&
                        kv.second == reverse(coefficients),
                graph.coeffs,
            ),
        )
        redundant_candidates =
            union(redundant_candidates_order, redundant_candidates_reverse)
    end

    # Rendundants nodes must both or neither be output nodes.
    is_output_node = node in graph.outputs
    redundant_nodes = collect(
        filter(
            x ->
                x in redundant_candidates &&
                    !xor(x in graph.outputs, is_output_node),
            redundant_candidates,
        ),
    )

    return redundant_nodes
end

# Merge redundant nodes into node.
function merge_redundant_nodes!(graph, node, redundant_nodes, cref, verbose)
    for key in redundant_nodes
        conditional_println("Merge node: $key with $node", verbose)
        replace_node!(graph, key, node, cref)
    end
end

"""
    compress_graph_redundant!(
        graph,
        cref = [];
        compress_lincomb = true,
        verbose = false,
    )

Removes from the graph redundant nodes, that is, nodes that repeat a computation
that is already present in the graph. Nodes corresponding to a linear
combination are removed only if the coefficients are the same and
`compress_lincomb` is set to `true`.
"""
function compress_graph_redundant!(
    graph,
    cref = [];
    compress_lincomb = true,
    verbose = false,
)
    # cref is not used as this function does not remove lincomb nodes.
    ismodified = true
    while ismodified
        ismodified = false
        for key in get_sorted_keys(graph) # No need to check input nodes.
            if key in get_sorted_keys(graph) # Check that node is not removed
                redundant_nodes =
                    find_redundant_nodes(graph, key, compress_lincomb)
                if !isempty(redundant_nodes)
                    merge_redundant_nodes!(
                        graph,
                        key,
                        redundant_nodes,
                        cref,
                        verbose,
                    )
                    ismodified = true
                end
            end
        end
    end
end

"""
    compress_graph_passthrough!(graph,cref=[];verbose=false);

Identifies lincombs lincomb of length one with coeff equal to one. The node has no effect. It redirect appropriately.
"""

function compress_graph_passthrough!(graph,
                                     cref=[];
                                     verbose = false)
    modified = true;
    while modified
        # all lincomb  keys
        lincomb_keys=collect(keys(graph.coeffs))
        lincomb_coeffs=map(x->graph.coeffs[x],lincomb_keys);

        # Find first occurance with one coeff and equal to 1
        i=findfirst(map(z-> z==[1], lincomb_coeffs))

        modified = false;
        if (!isnothing(i))
            modified = true
            node=lincomb_keys[i];
            parent=graph.parents[node][1];
            conditional_println(
                "Remove node: $node is replaced by $parent because of one mult lincomb", verbose
            )
            # node is computed using only one single scaling.
            GraphMatFun.replace_node!(graph,node,parent,cref);
        end
    end
    return graph
end

"""
    compress_graph!(graph,cref=[];verbose=false)

Searches for nodes which can be removed or reorganized without changing the
function it represents. Corresponding references in the `cref`-vector are
removed.
"""
function compress_graph!(graph, cref = []; verbose = false)
    compress_graph_output_cleaning!(graph, cref, verbose = verbose)
    compress_graph_zero_coeff!(graph, cref, verbose = verbose)
    compress_graph_trivial!(graph, cref, verbose = verbose)
    compress_graph_passthrough!(graph, cref, verbose = verbose)
    compress_graph_dangling!(graph, cref, verbose = verbose)
    return graph
end
