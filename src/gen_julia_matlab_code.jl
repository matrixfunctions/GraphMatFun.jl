## Code gen for julia and matlab

# Memory allocation (potentially with a new slot)
function get_mem!(mem,node,lang)
    if (haskey(mem,node))
        if (mem[node] isa Number)
            found_node=mem[node];
            if (lang==:matlab)
                return "memslots{$found_node}";
            else
                return "memslots[$found_node]";
            end
            delete!(mem,node);
        else
            return mem[node];
        end
    else # It's a new node. We need to find a free slot
        freeslot=-1;
        taken_slots=Set(values(mem))
        for k=1:1000 # Max number of memslots
            if (!(k in taken_slots))
                freeslot=k;
                break;
            end
        end
        if (freeslot == -1)
            error("No free slot found");
        end
        mem[node]=freeslot;
        if (lang==:matlab)
            return "memslots{$freeslot}";
        else
            return "memslots[$freeslot]";
        end
    end
end
function deallocate_mem!(mem,node,lang)
    found_node=-1;

    if (haskey(mem,node))
        if (mem[node] isa Number)
            found_node=mem[node];
            delete!(mem,node);
            if (lang==:matlab)
                return "memslots{$found_node}";
            else
                return "memslots[$found_node]";
            end
        else
            return mem[node];
        end
    else
        return "tmp_err";
    end
end

function write_operation_plain(file,T,graph,node,lang,op,
                               parent1,
                               parent2,
                               dealloc_list,
                               mem,mem_mode)


    nodemem=String(node);
    parent1mem=String(parent1);
    parent2mem=String(parent2);
    if op == :mult
        println(file,nodemem, "=",
                parent1mem, "*",
                parent2mem, ";")
    elseif op == :ldiv
        println(file, nodemem, "=",
                parent1mem, "\\",
                parent2mem, ";")
    elseif op == :lincomb
        for i=1:2
            print(file, "coeff$i=", real(graph.coeffs[node][i]))
            if T <: Complex
                t = imag(graph.coeffs[node][i])
                if t >= 0
                    print(file, " + ", t)
                else
                    print(file, " - ", abs(t))
                end
                print(file, "i")
            end
            println(file, ";")
        end
        if (graph.coeffs[node][1]==0)
            println(file, nodemem, "=",
                    "coeff2*",
                    parent2mem, ";")
        else
            println(file, nodemem, "=",
                    "coeff1*",
                    parent1mem, "+",
                    "coeff2*",
                    parent2mem, ";")
        end
    end
end

function write_operation_julia(file,T,graph,node,lang,op,
                               parent1,
                               parent2,
                               dealloc_list,
                               mem,mem_mode)




    parent1mem= get_mem!(mem,parent1,lang);
    parent2mem= get_mem!(mem,parent2,lang);


    nodemem="";
    if op == :mult

        nodemem=get_mem!(mem,node,lang);
        println(file,"mul!($nodemem,$parent1mem,$parent2mem);");
    elseif op == :ldiv
        nodemem=get_mem!(mem,node,lang);
        if (parent1 in dealloc_list)
            println(file,"# Improvable since $parent1=$parent1mem can be forgotten and this can be done in place");
        end
        println(file, nodemem, "=",
                parent1mem, "\\",
                parent2mem, ";")
    elseif op == :lincomb
        for i=1:2
            print(file, "coeff$i=", real(graph.coeffs[node][i]))
            if T <: Complex
                t = imag(graph.coeffs[node][i])
                if t >= 0
                    print(file, " + ", t)
                else
                    print(file, " - ", abs(t))
                end
                print(file, "i")
            end
            println(file, ";")
        end

        # This copy might be avoided if $parent2mem is not needed later
        if (parent1mem != "I" && parent2mem != "I")
            if (parent2 in dealloc_list)
                println(file,"# Smart lincomb without alloc ");
                nodemem=parent2mem;
                v=mem[parent2];
                mem[node]=v;
                delete!(mem,parent2);
                println(file,"BLAS.axpby!(coeff1,$parent1mem,coeff2,$nodemem);");
            elseif (parent1 in dealloc_list)
                println(file,"# Smart lincomb without alloc flip order ");
                nodemem=parent1mem;
                v=mem[parent1];
                mem[node]=v;
                delete!(mem,parent1);
                println(file,"BLAS.axpby!(coeff2,$parent2mem,coeff1,$nodemem);");

            else
                # We need to allocate a new slot
                nodemem=get_mem!(mem,node,lang);
                println(file,"$nodemem[:]=$parent2mem;");
                println(file,"BLAS.axpby!(coeff1,$parent1mem,coeff2,$nodemem);");
            end
        else
            nodemem=get_mem!(mem,node,lang);
            if (graph.coeffs[node][1]==0)
                println(file, nodemem, "[:]=",
                        "coeff2*",
                        parent2mem, ";")
            else
                println(file, nodemem, "[:]=",
                        "coeff1*",
                        parent1mem, "+",
                        "coeff2*",
                        parent2mem, ";")
            end
        end
    end

    return (nodemem,parent1mem,parent2mem);
end

function write_operation_matlab(file,T,node,nodename,op,parentname1,parentname2,dealloc_list)

    # Needs to be updated
    if op == :mult
        println(file,nodename, "=",
                parentname1, "*",
                parentname2, ";")
    elseif op == :ldiv
        println(file, nodename, "=",
                parentname1, "\\",
                parentname2, ";")
    elseif op == :lincomb
        for i=1:2
            print(file, "coeff$i=", real(graph.coeffs[node][i]))
            if T <: Complex
                t = imag(graph.coeffs[node][i])
                if t >= 0
                    print(file, " + ", t)
                else
                    print(file, " - ", abs(t))
                end
                print(file, "i")
            end
            println(file, ";")
        end
        if (graph.coeffs[node][1]==0)
            println(file, nodename, "=",
                    "coeff2*",
                    parentname2, ";")
        else
            println(file, nodename, "=",
                    "coeff1*",
                    parentname1, "+",
                    "coeff2*",
                    parentname2, ";")
        end
    end

end

#     mem_mode=:none, :prealloc, :clear
function gen_code(fname,graph;
                  priohelp=Dict{Symbol,Float64}(),
                  mem_mode=:prealloc,
                  lang=:julia,
                  funname="dummy")

    T=eltype(eltype(typeof(graph.coeffs.vals)));
    fname = abspath(fname)
    file = open(fname, "w+")



    cs="#"; # comment sign
    if (lang==:matlab)
        cs="%";
    end


    (order, can_be_deallocated, max_nof_nodes) =
        get_topo_order(graph; priohelp=priohelp);


    if (lang==:matlab)
        println(file,"function output=$funname(A);");
        println(file,"I=eye(size(A));");
    else
        println(file,"function $funname(A);");
        println(file,"i=1im;"); # Simplify complex numbers
    end




    if (mem_mode == :prealloc)

        println(file,"max_memslots=$max_nof_nodes + 1;");


        if (lang==:matlab)
            println(file,"$cs Allocate memory");
            println(file,"memslots=cell(max_memslots,1);");
            println(file,"for i=1:max_memslots");
            println(file,"   memslots{i}=zeros(size(A));");
            println(file,"end");
        else
            # Julia
            println(file,"T=promote_type(eltype(A),$T); $cs Make it work for many bigger types (matrices and scalars)");
            println(file,"memslots=Vector{Matrix{T}}(undef,max_memslots)");
            println(file,"n=size(A,1)");
            println(file,"for i=1:max_memslots");
            println(file,"    memslots[i]=Matrix{T}(undef,n,n);");
            println(file,"end");
        end
    end


    prealloc=max_nof_nodes+2;
    mem = Dict{Symbol,Union{Number,String}}();
    mem[:I]="I";
    mem[:A]="A";


    nof_written_outputs=0;
    for (i,node) in enumerate(order)
        op = graph.operations[node]
        parent1= graph.parents[node][1];
        parent2= graph.parents[node][2];

        nodemem= String(node);
        parent1mem= String(parent1);
        parent2mem= String(parent2);

        println(file,"$cs computing $node");


        if (lang==:matlab)
            if mem_mode==:prealloc # Memory alloc
                nodemem = get_mem!(mem,node,lang)
                parent1mem= get_mem!(mem,graph.parents[node][1],lang);
                parent2mem= get_mem!(mem,graph.parents[node][2],lang);
            end
            write_operation_matlab(file,T,node,nodemem,op,
                                   parent1mem,parent2mem,
                                   can_be_deallocated[i])
        else

            if (mem_mode == :none || mem_mode == :clear)
                write_operation_plain(file,T,graph,node,lang,op,
                                      parent1,
                                      parent2,
                                      can_be_deallocated[i],
                                      mem,mem_mode)

            else
                (nodemem,parent1mem,parent2mem)=
                write_operation_julia(file,T,graph,node,lang,op,
                                      parent1,
                                      parent2,
                                      can_be_deallocated[i],
                                      mem,mem_mode)
            end

        end


        if any(graph.outputs.==node)
            println(file, "output$nof_written_outputs", "=", nodemem, ";")
            println(file, "output = output$nof_written_outputs;");
            nof_written_outputs += 1
            if (lang==:julia)
                println(file,"return output");
            end

        end
        for n=can_be_deallocated[i]
            if ((n != :I) && (n!=:A))
                if (mem_mode==:prealloc)
                    #println("Deallocating $n");
                    nn=deallocate_mem!(mem,n,lang)
                    println(file,"$cs Effectively dealloc $n=$nn");
                elseif (mem_mode==:clear)
                    if (lang==:matlab)
                        println(file,"clear $n;");
                    else
                        println(file,"$n = [];");
                    end

                else
                    println(file,"$cs Could be cleared $n;");
                end

            end;
        end
        if (mem_mode==:prealloc)
            print(file,"$cs memlist=[");
            for (j,n) in enumerate(mem)
                if (j>1)
                    print(file,",");
                end
                print(file,"$n")
            end
            println(file,"]");;
        end


    end

    println(file,"end");
    close(file)
end
