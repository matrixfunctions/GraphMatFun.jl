# Main data structures of the computation graph,
# including sorting and manipulation functions

using LinearAlgebra

import Base.eltype
import Base.big
import Base.complex

export Compgraph
export get_topo_order
export add_mult!
export add_lincomb!
export add_ldiv!
export add_sum!
export add_output!
export del_output!
export clear_outputs!
export rename_node!
export del_node!

export get_all_cref
export extract_sums

export get_coeffs
export set_coeffs!

export get_sorted_keys

# Hash tables representing a computation graph
struct Compgraph{T}
    operations::Dict{Symbol,Symbol}
    parents::Dict{Symbol,Vector{Symbol}}
    coeffs::Dict{Symbol,Vector{T}}
    outputs::Vector{Symbol}
end

"""
    graph=Compgraph(T::Type=ComplexF64)

Creates an empty computation graph of with coefficients of type `T`.
"""
function Compgraph(T::Type = ComplexF64)
    return Compgraph(
        Dict{Symbol,Symbol}(),
        Dict{Symbol,Vector{Symbol}}(),
        Dict{Symbol,Vector{T}}(),
        Vector{Symbol}(),
    )
end
"""
    graph=Compgraph(T,orggraph::Compgraph)

Converts a graph such that the coefficients become type `T`.
Note that `T` can be a `Number` but can be other objects with defined
operations, or just `Any`.
"""
function Compgraph(T, orggraph::Compgraph)
    newcoeffs = Dict{Symbol,NTuple{<:Any,T}}()
    for node in keys(orggraph.coeffs)
        t = convert.(T,orggraph.coeffs[node])
        newcoeffs[node] = t
    end
    return Compgraph(
        orggraph.operations,
        orggraph.parents,
        newcoeffs,
        orggraph.outputs,
    )
end
# Convenience helpers
"""
    bgraph=big(graph::Compgraph{T})

Converts the coefficients of `graph` to type `big(T)` and returns `bgraph`,
which is of type `Compgraph{big(T)}`.

See also [`Compgraph(T,orggraph::Compgraph)`](@ref).
"""
function big(orggraph::Compgraph{T}) where {T}
    return Compgraph(big(T), orggraph)
end
"""
    cgraph=big(graph::Compgraph{T})

Converts the coefficients of `graph` to type `complex(T)` and returns `cgraph`,
which is of type `Compgraph{complex(T)}`.

See also [`Compgraph(T,orggraph::Compgraph)`](@ref).
"""
function complex(orggraph::Compgraph{T}) where {T}
    return Compgraph(complex(T), orggraph)
end

"""
    add_mult!(graph,node,p1,p2)

Adds a multiplication of node `p1` and `p2` to the
`graph`. The result is stored in node `node`.
"""
function add_mult!(graph, node, p1, p2)
    check_node_name_legality(graph, node)
    graph.operations[node] = :mult
    graph.parents[node] = [p1; p2];
    return nothing
end

"""
    add_lincomb!(graph,node,α1,p1,α2,p2)

The operation `α1*p1+α2*p2` is added to the graph.
The result is stored in node `node`.

See also [`add_sum!`](@ref).
"""
function add_lincomb!(graph, node, α1, p1, α2, p2)
    check_node_name_legality(graph, node)
    graph.operations[node] = :lincomb
    graph.parents[node] = [p1; p2]
    graph.coeffs[node] = [α1; α2]
    return nothing
end
"""
    add_lincomb!(graph,node,coeffs,nodes)

Adds a linear combination of the nodes (Vector or Tuple)
multiplied with coeffs (Vector or Tuple)
"""
function add_lincomb!(graph, node, coeffs, nodes)
    check_node_name_legality(graph, node)
    graph.operations[node] = :lincomb
    graph.parents[node] = collect(nodes)
    graph.coeffs[node] = collect(coeffs)
    return nothing
end


"""
    add_ldiv!(graph,node,p1,p2)

The operation `inv(p1)*p2` is added to the graph.
The result is stored in node `node`.
"""
function add_ldiv!(graph, node, p1, p2)
    check_node_name_legality(graph, node)
    graph.operations[node] = :ldiv
    graph.parents[node] = [p1; p2]
    return nothing
end

"""
    add_output!(graph,node)

Adds `node` to the bottom of the list of outputs of the graph.

See also [`eval_graph`](@ref), [`eval_jac`](@ref), and [`eval_runerr`](@ref).
"""
function add_output!(graph, node)
    # TODO: Add this?
    # for k = findall(graph.outputs .== node)
    #     popat!(graph.outputs,k)
    # end
    push!(graph.outputs, node)
    return nothing
end

"""
    del_output!(graph,node)

Removes the output `node` from the list of output nodes of the graph.
This function does not remove the node from the graph.
"""
function del_output!(graph, node)
    return deleteat!(graph.outputs, findall(x -> x == node, graph.outputs))
end

"""
    clear_outputs!(graph)

Clears the list of output nodes of the graph.
This function does not remove the node nodes from the graph.
"""
function clear_outputs!(graph)
    return empty!(graph.outputs)
end

"""
    rename_node!(graph,src,dest,cref=Vector())

This changes the name of the node `src` to `dest` and updates all references to
the node, including coefficient references in `cref` and `graph.outputs`.
"""
function rename_node!(graph, src, dest, cref = Vector())
    if src == dest
        return
    end

    # Check whether node exists in the graph, and whether is an input node.
    if src in keys(graph.parents)
        node_type = :full_node
    elseif src in Iterators.flatten(values(graph.parents))
        node_type = :input_node
    else
        error("Node $src not present in the graph.")
    end
    # Parents
    if node_type == :full_node
        graph.parents[dest] = graph.parents[src]
        delete!(graph.parents, src)
    end
    # Update also parent pointers
    for (key, value) in graph.parents
        newval = map(x-> x==src ? dest : x, value)
        graph.parents[key] = newval
    end
    if node_type == :input_node
        return
    end

    # Operations
    graph.operations[dest] = graph.operations[src]
    delete!(graph.operations, src)

    # Coeffs
    if (haskey(graph.coeffs, src))
        graph.coeffs[dest] = graph.coeffs[src]
        delete!(graph.coeffs, src)
    end

    # Output
    for (i, key) in enumerate(graph.outputs)
        if key == src
            graph.outputs[i] = dest
        end
    end

    # Crefs
    for (i, val) in enumerate(cref)
        if (val[1] == src)
            cref[i] = (dest, val[2])
        end
    end
    return nothing
end

"""
    cref=add_sum!(graph,node,c,nodelist,base_name=node)

Adds a linear combination of more than two nodes, given in `nodelist::Vector`,
with coefficients given in `c`. The `base_name` is temporary variables for the
summing. The sum is stored in `node`.

Returns `cref` list with references.
"""
function add_sum!(graph, node, c, nodelist, base_name = node)
    if size(nodelist, 1) == 1
        error("Summing one element not allowed.")
    elseif (size(nodelist, 1) == 2)
        # Direct call to lincomb if we sum two nodes
        add_lincomb!(graph, node, c[1], nodelist[1], c[2], nodelist[2])
        return [(node, 1), (node, 2)]
    end

    cref = Vector{Tuple{Symbol,Int}}()
    key = Symbol("$(base_name)2")
    add_lincomb!(graph, key, c[1], nodelist[1], c[2], nodelist[2])
    push!(cref, (key, 1))
    push!(cref, (key, 2))
    for k = 3:size(nodelist, 1)
        prev_key = key
        key = Symbol("$(base_name)$k")
        add_lincomb!(graph, key, 1, prev_key, c[k], nodelist[k])
        push!(cref, (key, 2))
    end
    # Rename the final sum as the target node
    rename_node!(graph, key, node, cref)
    return cref
end

"""
    del_node!(graph,node)

Deletes `node` from the graph, and all data associated with it, except for
`graph.outputs`.

See [`del_output!`](@ref).
"""
function del_node!(graph, node)
    delete!(graph.parents, node)
    delete!(graph.operations, node)
    delete!(graph.coeffs, node)
    # Leave the graph.outputs to user
    return nothing
end

function check_node_name_legality(graph, node)
    nodestr = String(node)
    re1 = r"^output[0-9]*$"
    readd = r".*\+.*"
    reldiv = r".*\\.*"
    remult = r".*\*.*"
    reequa = r".*=.*"
    if (
        (node == :graph_coeff_type) ||
        (node == :coeff1) ||
        (node == :coeff2) ||
        occursin(re1, nodestr) ||
        occursin(readd, nodestr) ||
        occursin(reldiv, nodestr) ||
        occursin(remult, nodestr) ||
        occursin(reequa, nodestr)
    )
        error(
            "Node name '$node' is a reserved keyword. " *
            "Hint: Names cannot have the form 'output<number>', " *
            "nor contain the symbols '+', '*', '\\', or '='.",
        )
    end
end

"""
    v=get_sorted_keys(graph)

Returns a list of all nodes, sorted.
"""
function get_sorted_keys(graph)
    return sort(collect(keys(graph.parents)))
end

"""
    v=get_all_cref(graph)

Returns a list with references to all coefficients in the graph.
The list is sorted.
"""
function get_all_cref(graph)
    kv = Vector{Tuple{Symbol,Int}}()
    for key = keys(graph.coeffs)
        for (i,_)=enumerate(graph.coeffs[key])
            push!(kv,(key,i));
        end
    end
    return kv
end

"""
    set_coeffs!(graph, x, cref=get_all_cref(graph))

Sets the coefficient values in the coefficients specified in `cref::Vector` to
the values in the vector in `x::Vector`.
"""
function set_coeffs!(graph, x, cref = get_all_cref(graph))
    if (cref isa Tuple) #Workaround for single coefficient
        cref = [cref]
    end
    if !(length(x) == length(cref))
        error(
            "Vector of coefficients and defined set of coefficients" *
            "do not have the same length.",
        )
    end
    for (idx, k) in enumerate(cref)
        node = k[1]
        parentnr = k[2]
        # Update only one element in the Tuple of coeffs.
        new_coeffs=collect(graph.coeffs[node]);
        new_coeffs[parentnr]=x[idx];
        graph.coeffs[node]=Tuple(new_coeffs);
    end
    return nothing
end



"""
    x=get_coeffs(graph, cref=get_all_cref(graph))

Gets the coefficient values for the coefficients specified in `cref::Vector`.
"""
function get_coeffs(graph, cref = get_all_cref(graph))
    if (cref isa Tuple) #Workaround for single coefficient
        cref = [cref]
    end
    T = eltype(graph)
    x = Vector{T}(undef, length(cref))
    for (idx, k) in enumerate(cref)
        x[idx] = graph.coeffs[k[1]][k[2]]
    end
    return x
end


## Misc helpers
"""
    T=eltpye(graph::Compgraph{T})

Returns the type `T` of a `graph`.
"""
function eltype(graph::Compgraph{T}) where {T}
    return T
end

#     children=get_children(graph,node)
# Returns the (direct) children of `node`.
function get_children(graph, node)
    c = Vector{Symbol}()
    for (child, parents) in graph.parents
        for s = 1:size(parents, 1)
            if (parents[s] == node)
                push!(c, child)
                break
            end
        end
    end
    return c
end

# Only for topo order
function nof_uncomputed_children(graph, node, vals)
    cv = get_children(graph, node)
    # Assumption: vals only contain keys of computed elements
    computed_nodes = keys(vals)
    setdiff!(cv, computed_nodes)
    return length(cv)
end

"""
    (order,can_be_deallocated,max_nodes)=get_topo_order(
        graph;
        priohelp = Dict{Symbol,Float64}(),
        free_mem_bonus = 1000,
        will_not_dealloc = [:I],
        input = :A,
    )

Computes a topological sort of `graph`, that is, an ordering of the nodes
following which the function the graph represents can be evaluated. The
`priohelp` kwarg can be used to obtain a different topological ordering by
changing the node priority. The code assumes that the nodes `:I` and `input` do
not need to be computed.

The code uses a heuristic to minimize pathwidth. It is based on a point system.
You can influence the computation order by providing a `priohelp`. If you want
node `:B4` to be computed earlier, you can set `priohelp[:B4]=-5000.0`. The
`free_mem_bonus` is used in the heuristic to prioritize the computation of nodes
which release other nodes. The vector `will_not_deallocate` influences the order
specifying nodes that will not be deallocated and therefore gets no
`free_mem_bonus`.

The return value `order` is a `Vector` of Symbols, and `can_be_deallocated` is a
`Vector{Vector{Symbol}}` where the element `i` specifies the `Symbols` that are
unused after step `i` in the ordering. The `max_nodes` is the pathwidth.
"""
function get_topo_order(
    graph;
    priohelp = Dict{Symbol,Float64}(),
    free_mem_bonus = 1000,
    will_not_dealloc = [:I],
    input = :A,
)
    # Assumed to be true for the computed nodes, and not exist for other nodes
    is_computed = Dict{Symbol,Bool}()
    is_computed[:I] = true
    is_computed[input] = true

    # TODO: Can this be optimized? If I is not needed? Saves memory...
    # To keep track of the maximum number of nodes alive at the same time
    is_still_needed = Dict{Symbol,Bool}()
    is_still_needed[:I] = true
    is_still_needed[input] = true
    max_nof_nodes = 2

    outputs = graph.outputs
    uncomputed = get_sorted_keys(graph)

    # To keep track of nodes that can be deallocated after computing <node>.
    can_be_deallocated = Vector{Vector{Symbol}}(undef, length(uncomputed))
    for n = 1:length(uncomputed)
        can_be_deallocated[n] = Vector{Symbol}(undef, 0)
    end

    comp_node_nr = 1
    computation_order = Vector{Symbol}()
    while (length(uncomputed) > 0)
        # make list of computable now
        now_computable = Dict{Symbol,Float64}()
        for node in uncomputed
            nof_parents = size(graph.parents[node], 1)
            if (all(map(x -> haskey(is_computed, x), graph.parents[node])))
                # Parents are computed, i.e., this node is currently computable.

                # Greedy heuristic pt 1: Check how many (uncomputed) children a
                # possible node has. Fewer = better since likely to "disappear
                # soon".
                nof_uc = nof_uncomputed_children(graph, node, is_computed)
                if (haskey(priohelp, node))
                    # Adjust for user input
                    nof_uc += priohelp[node]
                end
                now_computable[node] = nof_uc

                can_dealloc_parent = Vector{Bool}(undef, nof_parents)
                for i = 1:nof_parents
                    # Greedy heuristic pt 2: Adjust for deallocation
                    # possibilities of parents.
                    p = graph.parents[node][i]
                    can_dealloc_parent[i] = (
                        (nof_uncomputed_children(graph, p, is_computed) == 1) && !(any(outputs .== p))
                    )
                    # Update can_deallocate depending on `will_not_dealloc`

                    can_dealloc_parent[i] =
                        can_dealloc_parent[i] && !(p in will_not_dealloc)
                end

                if all(can_dealloc_parent) && !(allequal(graph.parents[node]))
                    # Prioritize if parent can be deallocated Set to
                    # -2*free_mem_bonus if two memory slots can deallocate after
                    # this.
                    now_computable[node] += -2 * free_mem_bonus
                elseif any(can_dealloc_parent)
                    # Set to -free_mem_bonus one memory slot can deallocate
                    # after this.
                    now_computable[node] += -free_mem_bonus
                end
            end
        end

        if !isempty(now_computable)
            # "compute" the node, i.e., update what is computed and uncomputed
            compute_node = argmin(now_computable) #Locally greedy selection
            nof_parents = size(graph.parents[compute_node], 1)

            is_computed[compute_node] = true
            setdiff!(uncomputed, [compute_node])
            push!(computation_order, compute_node)

            # Keep track of number of nodes alive.
            # Only parents to currently "computed" node are affected
            is_still_needed[compute_node] = true
            for i = 1:nof_parents
                parent = graph.parents[compute_node][i]
                if (nof_uncomputed_children(graph, parent, is_computed) == 0) &&
                   !(any(outputs .== parent)) # Nonoutput node, no children left
                    delete!(is_still_needed, parent)
                    push!(can_be_deallocated[comp_node_nr], parent)
                end
                if allequal(graph.parents[compute_node])
                    break
                end
            end
            comp_node_nr += 1
            nof_nodes = length(is_still_needed)
            if nof_nodes > max_nof_nodes
                max_nof_nodes = nof_nodes
            end
        elseif !isempty(uncomputed)
            error(
                "Graph is disconnected. Nodes ",
                uncomputed,
                " cannot be computed.",
            )
        end
    end
    return (computation_order, can_be_deallocated, max_nof_nodes)
end
