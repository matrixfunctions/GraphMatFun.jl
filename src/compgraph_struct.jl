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
    parents::Dict{Symbol,Tuple{Symbol,Symbol}}
    coeffs::Dict{Symbol,Tuple{T,T}}
    outputs::Vector{Symbol}
end


"""
    graph=Compgraph(T::Type=ComplexF64)

Creates an empty computation graph of with coefficients of type `T`.

    """
function Compgraph(T::Type=ComplexF64)
    return Compgraph(Dict{Symbol,Symbol}(), Dict{Symbol,Tuple{Symbol,Symbol}}(),
     Dict{Symbol,Tuple{T,T}}(), Vector{Symbol}()
     );
end
"""
    graph=Compgraph(T,orggraph::Compgraph)

Converts a graph such that the coefficients become type `T`.
Note that `T` can be a `Number` but can be other objects with defined
operations, or just `Any`.

    """
function Compgraph(T,orggraph::Compgraph)
    newcoeffs=Dict{Symbol,Tuple{T,T}}();
    for node in keys(orggraph.coeffs)
        t=(convert(T,orggraph.coeffs[node][1]),
           convert(T,orggraph.coeffs[node][2]));
        newcoeffs[node]=t;
    end
    return Compgraph(orggraph.operations,
                     orggraph.parents,
                     newcoeffs,
                     orggraph.outputs
                     );
end
# Convenience helpers
"""
    bgraph=big(graph::Compgraph{T})

Converts the coefficients of `graph` to type `big(T)` and returns `bgraph`,
which is of type `Compgraph{big(T)}`.
    """
function big(orggraph::Compgraph{T}) where {T}
    return Compgraph(big(T),orggraph)
end
"""
    cgraph=big(graph::Compgraph{T})

Converts the coefficients of `graph` to type `complex(T)` and returns `cgraph`,
which is of type `Compgraph{complex(T)}`.
    """
function complex(orggraph::Compgraph{T}) where {T}
    return Compgraph(complex(T),orggraph)
end


"""
    add_mult!(graph,node,p1,p2)

Adds a multiplication of node `p1` and `p2` to the
`graph`. The result is stored in node `node`.
    """
function add_mult!(graph,node,p1,p2)
    check_node_name_legality(graph,node)
    graph.operations[node]=:mult
    graph.parents[node]=(p1,p2);
end

"""
    add_lincomb!(graph,node,α1,p1,α2,p2)

The operation `α1*p1+α2*p2` is added to the graph.
The result is stored in node `node`.
    """
function add_lincomb!(graph,node,α1,p1,α2,p2)
    check_node_name_legality(graph,node)
    graph.operations[node]=:lincomb
    graph.parents[node]=(p1,p2);
    graph.coeffs[node]=(α1,α2);
end

"""
    add_ldiv!(graph,node,p1,p2)

The operation `inv(p1)*p2` is added to the graph.
The result is stored in node `node`.
    """
function add_ldiv!(graph,node,p1,p2)
    check_node_name_legality(graph,node)
    graph.operations[node]=:ldiv
    graph.parents[node]=(p1,p2);
end

"""
    add_output!(graph,node)

Adds an output `node` to the graph.
    """
function add_output!(graph,node)
    # TODO: Add this?
    # for k = findall(graph.outputs .== node)
    #     popat!(graph.outputs,k)
    # end
    push!(graph.outputs,node);
end

"""
    del_output!(graph,node)

Removes the output `node` from the list of output nodes of the graph.
This function does not remove the node from the graph.
"""
function del_output!(graph,node)
    deleteat!(graph.outputs,findall(x->x==node,graph.outputs));
end

"""
    clear_outputs!(graph)

Clears the list of output nodes of the graph.
This function does not remove the node nodes from the graph.
"""
function clear_outputs!(graph)
    empty!(graph.outputs)
end

"""
    rename_node!(graph,src,dest,cref=Vector())

This changes the name of the node `src` to `dest` and updating
all references to the node including coefficient
references in `cref`.
    """
function rename_node!(graph,src,dest,cref=Vector())
    if src == dest
        return
    end

    # Check whether node exists in the graph, and whether is an input node
    if src in keys(graph.parents)
        node_type=:full_node
    elseif src in getindex.(values(graph.parents),1) ||
        src in getindex.(values(graph.parents),2)
        node_type=:input_node
    else
        error("Node $src not present in the graph.")
    end
    # Parents
    if node_type==:full_node
        graph.parents[dest]=graph.parents[src];
        delete!(graph.parents,src);
    end
    # Update also parent pointers
    for (key,value)=graph.parents
        newval=[value[1];value[2]]
        for s=1:2
            if (value[s] == src)
                newval[s]=dest;
            end
        end
        graph.parents[key]=Tuple(newval);
    end
    if node_type==:input_node
        return
    end

    # Operations
    graph.operations[dest]=graph.operations[src];
    delete!(graph.operations,src);

    # Coeffs
    if (haskey(graph.coeffs,src))
        graph.coeffs[dest]=graph.coeffs[src];
        delete!(graph.coeffs,src);
    end

    # Output
    for (i,key)=enumerate(graph.outputs)
        if key==src
            graph.outputs[i]=dest;
        end
    end

    # Crefs
    for (i,val)=enumerate(cref)
        if (val[1] == src)
            cref[i]=(dest,val[2]);
        end
    end
end

"""
    cref=add_sum!(graph,node,c,nodelist,base_name=node)

Adds a linear combination of more than two nodes, given
in `nodelist::Vector`, with coefficients given in `c`.
The `base_name` is temporary variables for the summing.
The sum is stored in `node`.

Returns `cref` list with references.

"""
function add_sum!(graph,node,c,nodelist,base_name=node)
    if size(nodelist,1)==1
        error("Summing one element not allowed.")
    elseif (size(nodelist,1)==2)
        # Direct call to lincomb if we sum two nodes
        add_lincomb!(graph,node,c[1],nodelist[1],c[2],nodelist[2])
        return [(node,1),(node,2)]
    end

    cref=Vector{Tuple{Symbol,Int}}()
    key=Symbol("$(base_name)2")
    add_lincomb!(graph,key,
                 c[1], nodelist[1],
                 c[2], nodelist[2]);
    push!(cref,(key,1))
    push!(cref,(key,2))
    for k=3:size(nodelist,1)
        prev_key=key;
        key=Symbol("$(base_name)$k")
        add_lincomb!(graph,key,
                     1,prev_key,
                     c[k],nodelist[k]);
        push!(cref,(key,2));
    end
    # Rename the final sum as the target node
    rename_node!(graph,key,node,cref);
    return cref;
end

"""
    del_node!(graph,node)

Deletes a node from the graph, and all data associated with it, except
for `graph.outputs`.
    """
function del_node!(graph,node)
    delete!(graph.parents,node);
    delete!(graph.operations,node);
    delete!(graph.coeffs,node);
    # Leave the graph.outputs to user
end

# Add an artificial node defined as:
#   :ArtificialX = 1.0 * output[end] + 0.0 * node
# sets the output to :ArtificialX. This modification
# does not change the input-output relation of the graph.
# Returns the symbol of the new node.
function add_artificial!(graph,node)
    key=:None
    for k=1:100
        tkey=Symbol("Artificial$k");
        if (!haskey(graph.parents,tkey))
            key=tkey;
            break
        end
    end
    if (key==:None)
        error("No new Artificial key slot found");
    end
    if (node==key)
        error("Cannot add an artificial node it itself")
    end

    # Set to eps() to make it not disappear in plotting
    add_lincomb!(graph,key,1.0,graph.outputs[end],eps(),node);
    graph.outputs[end]=key;
    return key;
end

function check_node_name_legality(graph,node)
    nodestr = String(node)
    re1 = r"^output[0-9]*$";
    readd = r".*\+.*"
    reldiv = r".*\\.*"
    remult = r".*\*.*"
    reequa = r".*=.*"
    if (  (node==:graph_coeff_type) || (node==:coeff1) || (node==:coeff2) ||
          occursin(re1,nodestr) || occursin(readd,nodestr) ||
          occursin(reldiv,nodestr) || occursin(remult,nodestr) || occursin(reequa,nodestr) )
        error(string("Node name '", node, "' is a reserved keyword. ", #
               "Hint: Names of type 'output' + a number are not allowed, ",#
               "and neither are names containing the symbols '+', '*', '", "\\", "', or '='."))
    end
end

"""
    v=get_sorted_keys(graph)

Returns a list of all nodes, sorted.
    """
function get_sorted_keys(graph)
    return sort(collect(keys(graph.parents)));
end

"""
    v=get_all_cref(graph)

Returns a list with references to all coefficients in the graph. The list is sorted.
    """
function get_all_cref(graph)
    k = sort(collect(keys(graph.coeffs)))
    kv = Vector{Tuple{Symbol,Int}}(undef,2*length(k))
    for i =1:length(k)
        idx = 2*(i-1) + 1;
        kv[idx] = (k[i],1);
        kv[idx+1] = (k[i],2);
    end
    return kv;
end

"""
    v=extract_sums(graph)

::Vector{Tuple{Vector{Float64},Vector{Symbol},Vector{Symbol}}}`

Returns a representation of sums in `graph` which may potentially be merged.
The vector `sums` contains a tuple for each of these sums. The three entries of
the tuple are:

* a vector of `Float64` values that represent the coefficients of the summands;
* a vector of `Symbol`s that correspond to the summands; and
* a vector of intermediate `Symbol`s (i.e, nodes) that can be merged.

The first two vectors have the same number of entries, one for each element that
can be merged in the sum.
 """
function extract_sums(graph)
    coeff,nodes,merged,sums=find_mergeable_sums(graph,graph.outputs[1],[])
    return sums
end

# Return true if `node` has multiple parents. Such nodes cannot be freed.
function has_multiple_parents(graph,node)
    sum(map(x->any(x.==node),values(graph.parents))) > 1
end

function find_mergeable_sums(graph,node,processed,curr_coeff=1)
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
        return Float64[],Symbol[],Symbol[],[]
    else
        # Call function recursively on both parents.
        push!(processed,node)
        (parent1,parent2)=graph.parents[node]
        curr_lincomb=graph.operations[node] == :lincomb
        (coeff1,coeff2)=curr_lincomb ? graph.coeffs[node] : (1,1)
        pcoeffs1,pnodes1,pmerged1,sums1=find_mergeable_sums(graph,parent1,processed,coeff1)
        pcoeffs2,pnodes2,pmerged2,sums2=find_mergeable_sums(graph,parent2,processed,coeff2)
        if curr_lincomb
            # Grow the sum by adjoining terms coming from parents, if any.
            has_multiple_parents(graph,node) && (curr_coeff=1)
            new_coeffs=vcat(curr_coeff*pcoeffs1,curr_coeff*pcoeffs2);
            new_nodes=vcat(pnodes1,pnodes2)
            new_merged=vcat(pmerged1,pmerged2)
            # Lincomb parents are added to the mergeable nodes.
            # Non-lincomb parents are added to the sum and put, and their
            # coefficients are added to the vector of coefficients.
            if haskey(graph.operations,parent1) &&
                graph.operations[parent1] == :lincomb &&
                !has_multiple_parents(graph,parent1)
                new_merged=vcat(new_merged,parent1)
            else
                new_coeffs=vcat(new_coeffs,curr_coeff*coeff1);
                new_nodes=vcat(new_nodes,parent1)
            end
            if haskey(graph.operations,parent2) &&
                graph.operations[parent2] == :lincomb &&
                !has_multiple_parents(graph,parent2)
                new_merged=vcat(new_merged,parent2)
            else
                new_coeffs=vcat(new_coeffs,curr_coeff*coeff2);
                new_nodes=vcat(new_nodes,parent2)
            end
            if node in graph.outputs || has_multiple_parents(graph,node)
                return Float64[],Symbol[],Symbol[],
                vcat(sums1,sums2,(new_coeffs,new_nodes,vcat(new_merged,node)))
            else
                return new_coeffs,new_nodes,new_merged,vcat(sums1,sums2)
            end
        else
            # Current node is not a lincomb node.
            sums=vcat(sums1,sums2)
            # Add sum of lincomb parents, if any, to sums.
            sums=!isempty(pcoeffs1) ?
                vcat(sums,(pcoeffs1,pnodes1,vcat(pmerged1,parent1))) : sums
            sums=!isempty(pcoeffs2) ?
                vcat(sums,(pcoeffs2,pnodes2,vcat(pmerged2,parent2))) : sums
            return Float64[],Symbol[],Symbol[],sums
        end
    end
end

"""
    set_coeffs!(graph, x, cref=get_all_cref(graph))

Sets the coefficient values in the coefficients specified in `cref::Vector` to the values in the vector in `x::Vector`.
    """
function set_coeffs!(graph, x, cref=get_all_cref(graph))
            if (cref isa Tuple) #Workaround for single coefficient
        cref=[cref]
    end
    if !(length(x)==length(cref))
        error("Vector of coefficients not the same length as defined set of coefficients.")
    end
    for (idx, k)= enumerate(cref)
        node = k[1]
        parentnr = k[2]
        if parentnr == 1
            other = graph.coeffs[node][2];
            graph.coeffs[node]=(x[idx], other);
        elseif parentnr == 2
            other = graph.coeffs[node][1];
            graph.coeffs[node]=(other, x[idx]);
        else
            error("This is undefined: Unknown parent.")
        end
    end
    return nothing
end

"""
    x=get_coeffs(graph, cref=get_all_cref(graph))

Gets the coefficient values for the coefficients specified in `cref::Vector`.
    """
function get_coeffs(graph, cref=get_all_cref(graph))
    if (cref isa Tuple) #Workaround for single coefficient
        cref=[cref]
    end
    T=eltype(graph)
    x=Vector{T}(undef,length(cref))
    for (idx, k) = enumerate(cref)
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

"""
    children=get_children(graph,node)

Returns the (direct) children of `node`.

    """
function get_children(graph,node)
    c=Vector{Symbol}();
    for (child,parents) in graph.parents
        for s=1:size(parents,1)
            if (parents[s]==node)
                push!(c,child);
                break;
            end
        end
    end
    return c;
end

# Only for topo order
function nof_uncomputed_children(graph,node,vals)
    cv=get_children(graph,node);
    computed_nodes = keys(vals) # Assumption: vals only contain keys to computed elements
    setdiff!(cv,computed_nodes)
    return length(cv);
end

"""
    (order,can_be_deallocated,max_nodes)=get_topo_order(graph; priohelp=Dict{Symbol,Float64}(),free_mem_bonus=1000,will_not_deallocate=[:I],input=:A)

Computes a topological sort of `graph`, that is, an ordering of the
nodes following which the function the graph represents can be
evaluated. The `priohelp` kwarg can be used to obtain a different
topological ordering by changing the node priority. The code assumes
that the nodes `:I` and `input` do not need to be computed.

The code uses a heuristic to minimize pathwidth. It is based on a
point system. You can influence the computation order by providing
a `priohelp`. If you want node `:B4` to be computed earlier,
 you can set `priohelp[:B4]=-5000.0`.  The `free_mem_bonus`
is used in the heuristic to prioritize the computation
of nodes which release other nodes. The vector `will_not_deallocate`
influences the order specifying nodes that will not be deallocated
and therefore gets no `free_mem_bonus`.

The return value `order` is a `Vector` of Symbols, and
`can_be_deallocated` is a `Vector{Vector{Symbol}}` where
the element `i` specifies the `Symbols` that are unused
after step `i` in the ordering. The `max_nodes` is
the pathwidth.

    """
function get_topo_order(graph; priohelp=Dict{Symbol,Float64}(),
                        free_mem_bonus=1000,will_not_dealloc=[:I], input=:A)
    # Assumed to be true for the computed nodes, and not exist for other nodes
    is_computed=Dict{Symbol,Bool}();
    is_computed[:I]=true;
    is_computed[input]=true;

    # TODO: Can this be optimized? If I is not needed? Saves memory...
    # To keep track of the maximum number of nodes alive at the same time
    is_still_needed=Dict{Symbol,Bool}();
    is_still_needed[:I]=true;
    is_still_needed[input]=true;
    max_nof_nodes = 2;

    outputs = graph.outputs;
    uncomputed=get_sorted_keys(graph);

    # To keep track of which nodes can be deallocated after computation of <node>
    can_be_deallocated=Vector{Vector{Symbol}}(undef,length(uncomputed));
    for n =  1:length(uncomputed)
        can_be_deallocated[n]=Vector{Symbol}(undef,0);
    end

    comp_node_nr = 1;
    computation_order=Vector{Symbol}();
    while (length(uncomputed)>0)
        # make list of computable now
        now_computable=Dict{Symbol,Float64}();
        for node in uncomputed;
            nof_parents = size(graph.parents[node],1);
            parent1=graph.parents[node][1];
            parent2=graph.parents[node][2];
            if (all(map(x->haskey(is_computed,x),graph.parents[node])))
                # Parents are computed, i.e., this node is currently computable

                # Greedy heuristic pt 1: Check how many (uncomputed) children a
                # possible node has. Fewer = better since likely to "disappear soon"
                nof_uc = nof_uncomputed_children(graph,node,is_computed);
                if (haskey(priohelp,node))
                    # Adjust for user input
                    nof_uc += priohelp[node];
                end
                now_computable[node] = nof_uc;

                can_dealloc_parent=Vector{Bool}(undef,nof_parents);
                for i = 1:nof_parents
                    # Greedy heuristic pt 2: Adjust for deallocation possibilities of parents.
                    p=graph.parents[node][i];
                    can_dealloc_parent[i] =
                        ( (nof_uncomputed_children(graph,p,is_computed)==1)
                          && !(any(outputs.==p)) );
                    # Update can_deallocate depending on `will_not_dealloc`

                    can_dealloc_parent[i] = can_dealloc_parent[i] && !(p in will_not_dealloc)


                end


                if all(can_dealloc_parent) && !(parent1==parent2)
                    # Prioritize if parent can be deallocated
                    # Set to -2*free_mem_bonus if two memory slots can deallocate after this.
                    now_computable[node] += -2*free_mem_bonus;
                elseif any(can_dealloc_parent)
                    # Set to -free_mem_bonus one memory slot can deallocate after this.
                    now_computable[node] += -free_mem_bonus
                end
            end
        end

        if !isempty(now_computable)
            # "compute" the node, i.e., update what is computed and uncomputed
            compute_node=findmin(now_computable)[2]; #Locally greedy selection
            nof_parents = size(graph.parents[compute_node],1);

            is_computed[compute_node]=true;
            setdiff!(uncomputed,[compute_node]);
            push!(computation_order,compute_node);

            # Keep track of number of nodes alive.
            # Only parents to currently "computed" node are affected
            is_still_needed[compute_node] = true;
            for i = 1:nof_parents
                parent=graph.parents[compute_node][i];
                if (nof_uncomputed_children(graph,parent,is_computed)==0) &&
                  !(any(outputs.==parent)) #Parent has no children left and is not an output
                    delete!(is_still_needed,parent);
                    push!(can_be_deallocated[comp_node_nr],parent);
                end
                if graph.parents[compute_node][1] == graph.parents[compute_node][2]
                    break;
                end
            end
            comp_node_nr += 1;
            nof_nodes = length(is_still_needed);
            if nof_nodes > max_nof_nodes
                max_nof_nodes = nof_nodes;
            end
        elseif !isempty(uncomputed)
            error("Graph is disconnected. Nodes ", uncomputed, " cannot be computed.")
        end
    end
    return (computation_order, can_be_deallocated, max_nof_nodes)
end
