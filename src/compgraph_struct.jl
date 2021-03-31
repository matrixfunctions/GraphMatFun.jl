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
export rename_node!
export del_node!

export get_all_cref

export get_coeffs
export set_coeffs!

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

Converts a graph such that the coefficients become type `T`. Note that `T` can be a `Number` but can be other objects with defined operations, or just `Any`.

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
function big(orggraph::Compgraph{T}) where {T}
    return Compgraph(big(T),orggraph)
end
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

Adds an output `node` to the graph function.
    """
function add_output!(graph,node)
    # TODO: Add this?
    # for k = findall(graph.outputs .== node)
    #     popat!(graph.outputs,k)
    # end
    push!(graph.outputs,node);
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

    # Parents
    graph.parents[dest]=graph.parents[src];
    delete!(graph.parents,src);
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
    crefs=add_sum!(graph,node,c,nodelist,base_name=node)

Adds a linear combination of more than two nodes, given
in `nodelist::Vector`, with coefficients given in `c`.
The `base_name` is temporary variables for the summing.
The sum is stored in `node`.

Returns cref list with references.

"""
function add_sum!(graph,node,c,nodelist,base_name=node)
    if size(nodelist,1)==1
        error("Summing one element not allowed.")
    elseif (size(nodelist,1)==2)
        # Direct call to lincomb if we sum two nodes
        add_lincomb!(graph,node,c[1],nodelist[1],c[2],nodelist[2])
        return ((node,1),(node,2));
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
            error("This is undefine: Unknown parent.")
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
        for s=1:2
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
    (order,can_be_deallocated,max_nodes)=get_topo_order(graph; priohelp=Dict{Symbol,Float64}())

Computes a vector of all nodes sorted in a topological
way, i.e., an order it can be computed. The `priohelp`
kwarg can be used to obtain a different topological
ordering, by changing the node priority.

The code uses a heuristic to minimize pathwidth.

The return value `order` is a `Vector` of Symbols, and
`can_be_deallocated` is a `Vector{Vector{Symbol}}` where
the element `i` specifies the `Symbols` that are unused
after step `i` in the ordering. The `max_nodes` is
the pathwidth.

    """
function get_topo_order(graph; priohelp=Dict{Symbol,Float64}())
    # Assumed to be true for the computed nodes, and not exist for other nodes
    is_computed=Dict{Symbol,Bool}();
    is_computed[:I]=true;
    is_computed[:A]=true;

    # TODO: Can this be optimized? If I is not needed? Saves memory...
    # To keep track of the maximum number of nodes alive at the same time
    is_still_needed=Dict{Symbol,Bool}();
    is_still_needed[:I]=true;
    is_still_needed[:A]=true;
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
            parent1=graph.parents[node][1];
            parent2=graph.parents[node][2];
            if (haskey(is_computed,parent1) &&
                haskey(is_computed,parent2))
                # Parents are computed, i.e., this node is currently computable

                # Greedy heuristic pt 1: Check how many (uncomputed) children a
                # possible node has. Fewer = better since likely to "disappear soon"
                nof_uc = nof_uncomputed_children(graph,node,is_computed);
                if (haskey(priohelp,node))
                    # Adjust for user input
                    nof_uc += priohelp[node];
                end
                now_computable[node] = nof_uc;

                # Greedy heuristic pt 2: Adjust for deallocation possibilities of parents.
                can_dealloc_parent1 = ( (nof_uncomputed_children(graph,parent1,is_computed)==1)
                                        && !(any(outputs.==parent1)) );
                can_dealloc_parent2 = ( (nof_uncomputed_children(graph,parent2,is_computed)==1)
                                        && !(any(outputs.==parent2)) );
                # Set to "-Inf" if one parent can be deallocated, and
                if  can_dealloc_parent1 && can_dealloc_parent2 && !(parent1==parent2)
                    # Set to NaN if two memory slots can deallocate after this.
                    # Since, NaN is soreted before -Inf
                    now_computable[node] = NaN;
                elseif can_dealloc_parent1 || can_dealloc_parent2
                    # Set to -Inf one memory slot can deallocate after this.
                    # Since, NaN is soreted before -Inf
                    now_computable[node] = -Inf;
                end
            end
        end

        if !isempty(now_computable)
            # "compute" the node, i.e., update what is computed and uncomputed
            compute_node=findmin(now_computable)[2]; #Locally greedy selection
            is_computed[compute_node]=true;
            setdiff!(uncomputed,[compute_node]);
            push!(computation_order,compute_node);

            # Keep track of number of nodes alive.
            # Only parents to currently "computed" node are affected
            is_still_needed[compute_node] = true;
            for i = 1:2
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
