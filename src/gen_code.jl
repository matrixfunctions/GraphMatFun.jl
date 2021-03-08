include("gen_code_mem.jl");

export gen_code,LangJulia,LangMatlab,LangC

struct LangMatlab end;
struct LangJulia
    overwrite_input # Overwrite input
    exploit_uniformscaling  # I=UniformScaling. Should it be exploited?
end;
struct LangC end;


function LangJulia()
    return LangJulia(true,true);
end


# Code snippet handling
struct CodeSnippet
    code_lines::Vector{String}
    lang
end
function init_code(lang)
    return CodeSnippet(Vector{String}(undef,0),lang);
end
function push_code!(code,str)
    push!(code.code_lines,str);
end
function push_comment!(code,str)
    push!(code.code_lines,comment(code.lang,str));
end
function to_string(code)
    return join(code.code_lines,"\n");
end




# Every language needs:
# comment(::Lang,s)
# slotname(::Lang,i) #
# assign_coeff(::Lang,v,i)
# function_definition(::Lang,funname)
# function_init(lang::Lang,T,mem)
# init_mem(lang::Lang,max_nof_nodes)
# function_end(lang::Lang,graph,mem)
# execute_operation!(lang::Lang,T,graph,node,
#                    dealloc_list,   mem)


# Language specific comment operations
comment(::LangMatlab,s)="% $s";
comment(::LangJulia,s)="# $s";
comment(::LangC,s)="/* $s */";


slotname(::LangJulia,i)="memslots[$i]"
slotname(::LangMatlab,i)="memslots{$i}"

# For julia more parsing may be required TODO
assign_coeff(::LangJulia,v,i)=("coeff$i","coeff$i=$v");
assign_coeff(::LangMatlab,v,i)=("coeff$i","coeff$i=$(real(v)) + 1i*$(imag(v))");




### Julia

function function_definition(lang::LangJulia,funname)
    code=init_code(lang);
    push_code!(code,"using LinearAlgebra");
    push_code!(code,"function $funname(A)");
    return code
end
function function_init(lang::LangJulia,T,mem)
    code=init_code(lang);
    max_nodes=size(mem.slots,1);

    # Allocation
    push_code!(code,"max_memslots=$max_nodes;");
    push_code!(code,"T=promote_type(eltype(A),$T); "*comment(lang,"Make it work for many 'bigger' types (matrices and scalars)"))
    push_code!(code,"memslots=Vector{Matrix{T}}(undef,max_memslots)");
    push_code!(code,"n=size(A,1)");
    start_j=2;

    push_comment!(code,"The first slots are I and A");
    if (lang.overwrite_input)
        start_j=3
    end
    push_code!(code,"for  j=$start_j:max_memslots");

    push_code!(code,"    memslots[j]=Matrix{T}(undef,n,n);");
    push_code!(code,"end");
    # Initialize I
    I_slot_name=get_slot_name(mem,:I);
    if (!lang.exploit_uniformscaling)
        push_code!(code,"$I_slot_name=Matrix{T}(I,n,n)")
    else
        push_comment!(code,"Uniform scaling is exploited. No I matrix explicitly allocated");
    end

    # Initialize A
    A_slot_name=get_slot_name(mem,:A);
    if (lang.overwrite_input)
        # Overwrite input A
        push_code!(code,"$A_slot_name=A "*comment(lang,"overwrite A"))
    else
        # Otherwise make a copy
        push_code!(code,"$A_slot_name[:]=A ")
    end


    return code
end

function init_mem(lang::LangJulia,max_nof_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i));
    alloc_slot!(mem,1,:I);
    alloc_slot!(mem,2,:A);
    if (lang.exploit_uniformscaling)
        mem.special_names[:I]="I";
    end
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


    # Don't overwrite the list
    dealloc_list=deepcopy(dealloc_list);
    setdiff!(dealloc_list,keys(mem.special_names));
    parent1mem= get_slot_name(mem,parent1)
    parent2mem= get_slot_name(mem,parent2)

    code=init_code(lang);
    push_comment!(code,"Computing $node with operation: $op");
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
        coeff_vars=Vector{String}(undef,2);
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1);
        push_code!(code,coeff1_code);
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2);
        push_code!(code,coeff2_code);

        # Use economical memory slots
        if ((parent1 in dealloc_list) || (parent2 in dealloc_list))
            # Smart / inplace: parentX can be used to store newly computed value

            # No allocation needed

            recycle_parent=(parent1 in dealloc_list) ? parent1 : parent2

            push_comment!(code,"Smart lincomb recycle $recycle_parent");


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

        push_comment!(code,"Deallocating $n in slot $i");
        free!(mem,i);

    end

    return (code,nodemem)
end



### MATLAB
function function_definition(lang::LangMatlab,funname)
    code=init_code(lang);
    push_code!(code,"function output=$funname(A)");
    return code
end

function function_init(lang::LangMatlab,T,mem)
    code=init_code(lang);
    push_code!(code,"I=eye(n,n);");
    return code;
end

function init_mem(lang::LangMatlab,max_nof_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i));
    return mem;
end
function function_end(lang::LangMatlab,graph,mem)
    code=init_code(lang);
    push_code!(code,"output=$(graph.outputs[end])");
    push_code!(code,"end")
    return code
end


function execute_operation!(lang::LangMatlab,
                            T,graph,node,
                            dealloc_list,
                            mem)
    op = graph.operations[node]
    parent1=graph.parents[node][1]
    parent2=graph.parents[node][2]


    code=init_code(lang);
    push_comment!(code,"Computing $node with operation: $op");
    # Needs to be updated
    if op == :mult
        push_code!(code,"$node=$parent1 * $parent2;");
    elseif op == :ldiv
        push_code!(code,"$node=$parent1 \\ $parent2;");
    elseif op == :lincomb
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1);
        push_code!(code,"$coeff1_code;");
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2);
        push_code!(code,"$coeff2_code;");
        push_code!(code,"$node= $coeff1*$parent1+$coeff2*parent2;")
    end
    return (code,"$node");
end







## Main code generation function
"""
    gen_code(fname,graph; priohelp=Dict{Symbol,Float64}(),
             lang=LangJulia(),funname="dummy")

Generates the code for the `graph` in the language
 specified in `lang` and writes it into the file
`fname`. The string `funname` is the function name.
Topological order of the nodes is comptued using
`get_topo_order` and `priohelp` can be used to
influence the order.

Currently supported language: `LangJulia`, `LangMatlab`, `LangC`.

"""
function gen_code(fname,graph;
                  priohelp=Dict{Symbol,Float64}(),
                  lang=LangJulia(),
                  funname="dummy")

    T=eltype(eltype(typeof(graph.coeffs.vals)));

    if (fname isa String)
        fname = abspath(fname)
        file = open(fname, "w+")
    else
        # Lazy: Print out to stdout if no filename
        file=Base.stdout;
    end

    (order, can_be_deallocated, max_nof_slots) =
        get_topo_order(graph; priohelp=priohelp);

    # max_nof_slots is the path width which gives
    # a bound on the number memory slots needed.


    println(file,to_string(function_definition(lang,funname)));

    # We do a double sweep in order to determine exactly
    # how many memory slots are needed. The first sweep
    # we carry out all operations but store only the
    # maximum number of memory slots needed. The
    # second sweep generates the code.


    mem=init_mem(lang,max_nof_slots+2)

    # Sweep 1: Determine exactly the number of slots needed
    nof_slots=0;
    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                     T,graph,node,
                                     can_be_deallocated[i],
                                     mem)

        # How many slots needed to reach this point
        if (!isnothing(findlast(mem.slots .!= :Free)))
            nof_slots=max(nof_slots,findlast(mem.slots .!= :Free))
        end

    end


    # Sweep 2:
    mem=init_mem(lang,nof_slots)
    println(file,to_string(function_init(lang,T,mem)));


    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                     T,graph,node,
                                     can_be_deallocated[i],
                                     mem)
        println(file,to_string(exec_code))
    end
    println(file,to_string(function_end(lang,graph,mem)));

    if (fname isa String)
        close(file)
    end


end
