export LangJulia

# Data structure for the language.
struct LangJulia
    overwrite_input # Overwrite input
    exploit_uniformscaling  # I=UniformScaling. Should it be exploited?
end
function LangJulia()
    return LangJulia(true,true)
end



# Language specific operations.
comment(::LangJulia,s)="# $s"

slotname(::LangJulia,i)="memslots[$i]"

assign_coeff(lang::LangJulia,v,i)=
    v==1 ?
    ("value_one",
     comment(lang,"Saving scalar multiplications using ValueOne().")) :
         ("coeff$i","coeff$i=$v")

# Code generation.
function push_code_matfun_axpy!(code)

    # Add code for matfun_axpby! functions.
    matfun_axpby_functions="
struct ValueOne; end
ValueOne()

# Compute X <- a X + b Y.
function matfun_axpby!(X,a,b,Y::UniformScaling)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:n
        X[i,i]+=(b isa ValueOne) ? 1 : b
    end
end
function matfun_axpby!(X,a,b,Y)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:m
        @inbounds for j=1:n
            if (b isa ValueOne)
                X[i,j]+=Y[i,j]
            else
                X[i,j]+=b*Y[i,j]
            end
        end
    end
end"
    push_code_verbatim_string!(code,matfun_axpby_functions)
end

function function_definition(lang::LangJulia,graph,T,funname)
    code=init_code(lang)
    push_code!(code,"using LinearAlgebra",ind_lvl=0)
    if any(values(graph.operations) .== :lincomb)
        push_code_matfun_axpy!(code)
    end
    push_code!(code,"function $funname(A)",ind_lvl=0)
    return code
end

function function_init(lang::LangJulia,T,mem,graph)
    code=init_code(lang)
    max_nodes=size(mem.slots,1)

    # Allocation
    push_code!(code,"max_memslots=$max_nodes")
    push_code!(code,"T=promote_type(eltype(A),$T) "*
        comment(lang,"Make it work for many 'bigger' types (matrices and scalars)"))
    push_code!(code,"memslots=Vector{Matrix{T}}(undef,max_memslots)")
    push_code!(code,"n=size(A,1)")
    start_j=2

    push_comment!(code,"The first slots are I and A")
    if (lang.overwrite_input)
        start_j=3
    end
    push_code!(code,"for  j=$start_j:max_memslots")

    push_code!(code,"memslots[j]=Matrix{T}(undef,n,n)",ind_lvl=2)
    push_code!(code,"end")

    # If needed, initialize variable of type ValueOne for axpby with a=1 or b=1.
    if (any(map(x->any(x .== 1),values(graph.coeffs))))
        push_code!(code,"value_one=ValueOne()")
    end

    # Initialize I
    I_slot_name=get_slot_name(mem,:I)
    if (!lang.exploit_uniformscaling)
        push_code!(code,"$I_slot_name=Matrix{T}(I,n,n)")
    else
        push_comment!(code,"Uniform scaling is exploited. No I matrix explicitly allocated")
    end

    # Initialize A.
    A_slot_name=get_slot_name(mem,:A)
    if (lang.overwrite_input)
        # Overwrite input A
        push_code!(code,"$A_slot_name=A "*comment(lang,"overwrite A"))
    else
        # Otherwise make a copy
        push_code!(code,"copy!($A_slot_name,A)")
    end

    return code
end

function init_mem(lang::LangJulia,max_nof_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i))
    alloc_slot!(mem,1,:I)
    alloc_slot!(mem,2,:A)
    if (lang.exploit_uniformscaling)
        mem.special_names[:I]="I"
    end
    return mem
end

function function_end(lang::LangJulia,graph,mem)
    code=init_code(lang)
    retval_node=graph.outputs[end]
    retval=get_slot_name(mem,retval_node)

    push_code!(code,"return $retval "*comment(lang,"Returning $retval_node"))
    push_code!(code,"end",ind_lvl=0)
    return code
end

# Adding of a scaled identity matrix in julia code.
function execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
    push_comment!(code,"Add lincomb with identity.")

    if (nodemem != non_I_parent_mem)
        push_code!(code,"copy!($(nodemem),$(non_I_parent_mem))")
    else
        push_comment!(code,"No copy necessary. Inline identity multiple add")
    end

    # push_code!(code,"$(nodemem) .*= $non_I_parent_coeff")
    # push_code!(code,"D=view($(nodemem), diagind($(nodemem), 0))")
    # push_code!(code,"D .+= $I_parent_coeff")
    push_code!(code,"matfun_axpby!($(nodemem),$non_I_parent_coeff,$I_parent_coeff,I)")
end


function execute_operation!(lang::LangJulia,
                            T,graph,node,
                            dealloc_list,
                            mem)

    op = graph.operations[node]
    parent1=graph.parents[node][1]
    parent2=graph.parents[node][2]

    # Don't overwrite the list
    dealloc_list=deepcopy(dealloc_list)
    setdiff!(dealloc_list,keys(mem.special_names))
    parent1mem=get_slot_name(mem,parent1)
    parent2mem=get_slot_name(mem,parent2)

    code=init_code(lang)
    push_comment!(code,"Computing $node with operation: $op")
    if op == :mult
        # Multiplication has no inline
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,node)

        push_code!(code,"mul!($nodemem,$parent1mem,$parent2mem)")

    elseif op == :ldiv
        # Left division
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,node)
        push_code!(code,"$nodemem=$parent1mem\\$parent2mem")

    elseif op == :lincomb
        coeff_vars=Vector{String}(undef,2)
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1)
        push_code!(code,coeff1_code)
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2)
        push_code!(code,coeff2_code)

        # Use economical memory slots
        if ((parent1 in dealloc_list) || (parent2 in dealloc_list))
            # Smart / inplace: parentX can be used to store newly computed value

            # No allocation needed

            recycle_parent=(parent1 in dealloc_list) ? parent1 : parent2

            push_comment!(code,"Smart lincomb recycle $recycle_parent")

            nodemem=get_slot_name(mem,recycle_parent)

            if (lang.exploit_uniformscaling &&  (parent1 == :I || parent2 == :I))
                # BLAS does not work with unform scaling use inplace instead
                if (parent1 == :I)
                    non_I_parent_mem=parent2mem
                    non_I_parent_coeff=coeff2
                    I_parent_coeff=coeff1
                elseif (parent2 == :I)
                    non_I_parent_mem=parent1mem
                    non_I_parent_coeff=coeff1
                    I_parent_coeff=coeff2
                end
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
            else
                if (recycle_parent == parent1)
                    # push_code!(code,"BLAS.axpby!($coeff2,$parent2mem,$coeff1,$nodemem);")
                    push_code!(code,"matfun_axpby!($nodemem,$coeff1,$coeff2,$parent2mem)")
                else
                    # push_code!(code,"BLAS.axpby!($coeff1,$parent1mem,$coeff2,$nodemem);")
                    push_code!(code,"matfun_axpby!($nodemem,$coeff2,$coeff1,$parent1mem)")
                end
            end

            # Avoid deallocated
            setdiff!(dealloc_list,[recycle_parent])

            # Set the memory pointer
            j=get_slot_number(mem,recycle_parent)
            set_slot_number!(mem,j,node)
        else
            # Default behavior:  Allocate new slot
            (nodemem_i,nodemem)=get_free_slot(mem)

            alloc_slot!(mem,nodemem_i,node)

            if (parent1 == :I)
                non_I_parent_mem=parent2mem
                non_I_parent_coeff=coeff2
                I_parent_coeff=coeff1
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
            elseif (parent2 == :I)
                non_I_parent_mem=parent1mem
                non_I_parent_coeff=coeff1
                I_parent_coeff=coeff2
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
            else
                # push_code!(code,
                #            "$(nodemem)[:]=$coeff1*$parent1mem +
                # $coeff2*$parent2mem")
                push_code!(code,"copy!($(nodemem),$parent1mem)") # Arbitrary choice
                push_code!(code,"matfun_axpby!($(nodemem),$coeff1,$coeff2,$parent2mem)")
            end

        end

    else
        error("Unknown operation")
    end

    # Deallocate
    for n=dealloc_list
        i=get_slot_number(mem,n)

        push_comment!(code,"Deallocating $n in slot $i")
        free!(mem,i)
    end

    return (code,nodemem)
end
