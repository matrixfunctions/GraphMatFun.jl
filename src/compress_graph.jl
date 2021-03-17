# Functions to automatically remove
# nodes from graphs

export compress_graph_dangling!
export compress_graph_zero_coeff!
export compress_graph_output_cleaning!
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
end
