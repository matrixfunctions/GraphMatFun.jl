using Printf
include("gen_code_mem.jl");

export gen_code,LangJulia,LangMatlab,LangC

struct LangMatlab end;
struct LangJulia end;
struct LangC end;



# Language specific operations
comment(::LangMatlab,s)="% $s";
comment(::LangJulia,s)="# $s";
comment(::LangC,s)="/* $s */";


slotname(::LangJulia,i)="memslots[$i]"
slotname(::LangMatlab,i)="memslots{$i}"

# For julia more parsing may be required
assign_coeff(::LangJulia,v,i)=("coeff$i","coeff$i=$v");



function init_code(lang)
    return Vector{String}(undef,0);
end

function push_code!(code,str)
    push!(code,str);
end
function push_comment!(code,lang,str)
    push!(code,comment(lang,str));
end

function function_definition(::LangMatlab,funname)
    code=init_code(lang);
    push_code!(code,"function output=$funname(A)");
    return code
end
function function_definition(lang::LangJulia,funname)
    code=init_code(lang);
    push_code!(code,"function $funname(A)");
    return code
end
function function_init(lang::LangJulia,T,max_nodes)
    code=init_code(lang);
    push_code!(code,"i=1im");
    push_code!(code,"max_memslots=$max_nodes;");
    push_code!(code,"T=promote_type(eltype(A),$T); "*comment(lang,"Make it work for many 'bigger' types (matrices and scalars)"))


    push_code!(code,"memslots=Vector{Matrix{T}}(undef,max_memslots)");
    push_code!(code,"n=size(A,1)");
    push_code!(code,"for i=1:max_memslots");
    push_code!(code,"    memslots[i]=Matrix{T}(undef,n,n);");
    push_code!(code,"end");
    return code
end

function init_mem(lang::LangJulia,max_nof_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i));
    alloc_slot!(mem,1,:I);
    alloc_slot!(mem,2,:A);
    return mem;
end
function function_end(lang::LangJulia,graph,mem)
    code=init_code(lang);
    retval_node=graph.outputs[end];
    retval=get_slot_name(mem,retval_node);

    push_code!(code,"return $retval; "*comment(lang,"Returning $retval_node"))
    push_code!(code,"end");
    return code
end


function execute_operation!(lang::LangJulia,
                            T,graph,node,
                            dealloc_list,
                            mem)


    op = graph.operations[node]
    parent1=graph.parents[node][1]
    parent2=graph.parents[node][2]


    dealloc_list=deepcopy(dealloc_list); # Don't overwrite
    parent1mem= get_slot_name(mem,parent1)
    parent2mem= get_slot_name(mem,parent2)

    code=init_code(lang);
    push_comment!(code,lang,"Computing $node with operation: $op");
    if op == :mult
        # Multiplication has no inline
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,node);


        push_code!(code,"mul!($nodemem,$parent1mem,$parent2mem);");

    elseif op == :ldiv
        # Left division
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,node);
        push_code!(code,"$nodemem=$parent1mem\\$parent2mem")

    elseif op == :lincomb
        # Use economical memory slots
        coeff_vars=Vector{String}(undef,2);
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1);
        push_code!(code,coeff1_code);
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2);
        push_code!(code,coeff2_code);

        if ((parent1 in dealloc_list) || (parent2 in dealloc_list))
            # Smart / inplace: parentX can be used to store newly computed value

            # No allocation needed

            recycle_parent=(parent1 in dealloc_list) ? parent1 : parent2

            push_comment!(code,lang,"Smart lincomb recycle $recycle_parent");


            nodemem=get_slot_name(mem,recycle_parent)

            if (recycle_parent == parent1)
                push_code!(code,"BLAS.axpby!($coeff2,$parent2mem,$coeff1,$nodemem);");
            else
                push_code!(code,"BLAS.axpby!($coeff1,$parent1mem,$coeff2,$nodemem);");
            end

            # Avoid deallocated
            setdiff!(dealloc_list,[recycle_parent]);

            # Set the memory pointer
            j=get_slot_number(mem,recycle_parent);
            set_slot_number!(mem,j,node);

        else
            # Default behavior:  Allocate new slot
            (nodemem_i,nodemem)=get_free_slot(mem)

            alloc_slot!(mem,nodemem_i,node);
            push_code!(code,
                       "$(nodemem)[:]=$coeff1*$parent1mem + $coeff2*$parent2mem")
        end

    else
        error("Unknown operation");
    end


    # Deallocate
    for n=dealloc_list
        i=get_slot_number(mem,n);

        push_comment!(code,lang,"Deallocating $n in slot $i");
        free!(mem,i);

    end

    return (code,nodemem)
end





function gen_code(fname,graph;
                  priohelp=Dict{Symbol,Float64}(),
                  mem_mode=:prealloc,
                  lang=LangJulia(),
                  funname="dummy")

    T=eltype(eltype(typeof(graph.coeffs.vals)));

    file=Base.stdout;
    (order, can_be_deallocated, max_nof_nodes) =
        get_topo_order(graph; priohelp=priohelp);



    println(file,join(function_definition(lang,funname),"\n"));

    mem=init_mem(lang,max_nof_nodes+2)
    orgmem=deepcopy(mem);


    # Sweep 1: Determine the effective path width
    path_width=0;
    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                     T,graph,node,
                                     can_be_deallocated[i],
                                     mem)

        path_width=max(path_width,findlast(mem.slots .!= :Free))
    end


    # Sweep 2:
    println(file,comment(lang,"Path width $path_width"));

    println(file,join(function_init(lang,T,path_width),"\n"));
    mem=init_mem(lang,path_width)

    println(file,get_slot_name(mem,:I),"[:]=Matrix(I,n,n)")

    # Overwrite input A
    println(file,get_slot_name(mem,:A),"=A")

    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                     T,graph,node,
                                     can_be_deallocated[i],
                                     mem)
        println(file,join(exec_code,"\n"))
    end
    println(file,join(function_end(lang,graph,mem),"\n"));

end
