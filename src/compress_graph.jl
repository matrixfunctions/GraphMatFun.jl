# Functions to automatically remove
# nodes from graphs

export compress_graph_dangling!
export compress_graph_zero_coeff!
export compress_graph_output_cleaning!
export comprrss_graph_trivial_nodes!
export compress_graph!;


function delete_crefs!(cref,n)

    II=[]; # II will be list of nodes not to remove
    for (i,c)=enumerate(cref)
        if (c[1] != n)
            push!(II,i)
        end
    end
    cref_org=copy(cref);
    empty!(cref);
    append!(cref,cref_org[II])
end

"""
    compress_graph_dangling!(graph,cref=[])

Removes dangling nodes in the graph. If a node is not used anywhere
and it is not in the output, it can safely be removed.

    """
function compress_graph_dangling!(graph,cref=[])
    # Remove dangling nodes
    modified=true
    while modified
        modified=false;
        for (key,parents) in graph.parents
            nof_children=0;
            for (pkey,pparents) in graph.parents
                if (pparents[1]==key ||
                    pparents[2]==key)
                    nof_children += 1;
                end
            end

            if (nof_children==0 && !(key in graph.outputs))
                println("Delete dangling:",key);
                delete!(graph.parents,key)
                delete!(graph.operations,key)
                if (haskey(graph.coeffs,key))
                    delete!(graph.coeffs,key)
                end
                delete_crefs!(cref,key)
                modified=true;
                break
            end
        end
    end
end


"""
    compress_graph_zero_coeff!(graph,cref=[];droptol=0)

Searches for linear combinations with zero coeff value and tries to compress
by redirecting node references.

    """
function compress_graph_zero_coeff!(graph,cref=[];droptol=0)

    modified=true;
    while (modified)
        modified=false;
        for (key,parents) in graph.parents
            big_break=false;
            if (graph.operations[key] == :lincomb &&
                !(key in graph.outputs)) # don't remove outputs
                v=graph.coeffs[key]
                for s=1:2
                    other_s=3-s;
                    if (abs(v[s]) <= droptol)

                        println("Zero coeff can be removed: $key because of coeff $s is "*string(v[s]));

                        # Redirect all lincombs we can find
                        for (pkey,op) in graph.operations
                            pp=[graph.parents[pkey]...]
                            if (op==:lincomb &&
                                (pp[1]==key ||  pp[2]==key))
                                # It can be merged

                                println("Delete connection:",key,"->",pkey)
                                # z: the parent connected to the node we want to remove
                                z=1;
                                if (pp[2]==key)
                                    z=2;
                                end

                                cc=[graph.coeffs[pkey]...]
                                cc[z] *= cc[z]*v[other_s];

                                pp[z]=graph.parents[key][other_s];
                                # Overwrite
                                add_lincomb!(graph,pkey,
                                             cc[1],pp[1],
                                             cc[2],pp[2]);

                                modified=true;
                                big_break=true
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
    compress_graph_output_cleaning!(graph,cref=[])

Checks if the output is computed from the linear combination
that can be compressed.

    """
function compress_graph_output_cleaning!(graph,cref=[])

    ismodified=false;
    while ismodified
        ismodified=false;
        for (i,n)=enumerate(graph.outputs)
            if (graph.operations[n]==:lincomb)

                for s=1:2
                    other_s=3-s;
                    if (graph.coeffs[n][s]==1.0
                        &&
                        graph.coeffs[n][other_s]==0.0)
                        # Redirect since it is just a copy of
                        # another node
                        graph.outputs[i]=graph.parents[n][s];
                        ismodified=true;
                    end
                end
            end

        end
    end

end

### Remove trivial nodes:
# 1) :mult with identity parent;
# 2) :ldiv with identity as first (left) parent.

# Trivial node with identity on the left (:mult or :ldiv).
is_trivial_left(graph,node)=(graph.parents[node][1] == :I &&
    (graph.operations[node] == :mult || graph.operations[node] == :ldiv))
# Trivial node with identity on the right (:mult).
is_trivial_right(graph,node)=(graph.parents[node][2] == :I &&
    graph.operations[node] == :mult)
# Update children of oldparent to point to newparent.
function updated_children_deleted_node!(graph,oldparent,newparent)
    for child in get_children(graph,oldparent)
        graph.parents[child] = map(x->(x==oldparent ? newparent : x),
                                   graph.parents[child])
    end
end
"""
    has_identity_lincomb(graph)

Checks whether the graph has a node which is a linear combination of identity
matrices.

    """
function has_identity_lincomb(graph)
    has_identity_lincomb=false
    for (key,parents) in graph.parents
        if graph.operations[key]==:lincomb && parents[1]==:I && parents[2]==:I
            has_identity_lincomb=true
            break
        end
    end
    return has_identity_lincomb
end
"""
    has_trivial_node(graph)

Checks whether the graph has trivial nodes, that is, multiplications by the
identity or linear systems whose coefficient is the identity matrix.

    """
function has_trivial_nodes(graph)
    has_trivial_nodes=false
    for (key,parents) in graph.parents
        if is_trivial_left(graph,key) || is_trivial_right(graph,key)
            has_trivial_nodes=true
            break
        end
    end
    return has_trivial_nodes
end
"""
    compress_graph_trivial!(graph,cref=[])

Removes from the graph trivial nodes, that is, multiplications by the identity
or linear systems whose coefficient is the identity matrix.

    """
function compress_graph_trivial!(graph)
    ismodified=true
    while ismodified
        ismodified=false
        for (key,parents) in graph.parents
            delete_node=false
            if is_trivial_left(graph,key)
                newparent=graph.parents[key][2]
                delete_node=true
            elseif is_trivial_right(graph,key)
                newparent=graph.parents[key][1]
                delete_node=true
            end
            if (delete_node)
                println("Replace node ",key," by ",newparent);
                # Make newparent an output node if current node was.
                if (key in graph.outputs)
                    deleteat!(graph.outputs,graph.outputs .== key)
                    add_output!(graph,newparent)
                end
                update_children_deleted_node!(graph,key,newparent)
                del_node!(graph,key)
                ismodified=true
            end
        end
    end
end
"""
    compress_graph!(graph,cref=[])

Searches for nodes which can be removed or reorganized without changing
the function it represents. Corresponding references in the `cref`-vector
are removed.

    """
function compress_graph!(graph,cref=[])
    compress_graph_output_cleaning!(graph,cref)
    compress_graph_zero_coeff!(graph,cref)
    compress_graph_dangling!(graph,cref)
    compress_graph_trivial!(graph)
end
