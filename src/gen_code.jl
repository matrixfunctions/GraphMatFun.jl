include("gen_code_mem.jl");

export gen_code,LangJulia,LangMatlab,LangC

struct LangMatlab end;
struct LangJulia
    overwrite_input # Overwrite input
    exploit_uniformscaling  # I=UniformScaling. Should it be exploited?
end;
struct LangC_MKL end
struct LangC_OpenBLAS end
LangC=Union{LangC_MKL,LangC_OpenBLAS}

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
function push_code!(code,str;ind_lvl=1,ind_str="    ")
    # Indent only non-empty lines.
    indentation=repeat(ind_str,ind_lvl)
    indented_string=isempty(str) ? "" : indentation*str
    push!(code.code_lines,indented_string);
end
function push_comment!(code,str;ind_lvl=1,ind_str="    ")
    # Convert empty comments to empty lines.
    push_code!(code, isempty(str) ? "" : comment(code.lang,str),
               ind_lvl=ind_lvl,ind_str=ind_str)
end
function to_string(code)
    return join(code.code_lines,"\n");
end

# Every language needs:
# comment(::Lang,s)
# slotname(::Lang,i) #
# assign_coeff(::Lang,v,i)
# function_definition(::Lang,T,funname)
# function_init(lang::Lang,T,mem,graph)
# init_mem(lang::Lang,max_nof_nodes)
# function_end(lang::Lang,graph,mem)
# execute_operation!(lang::Lang,T,graph,node,
#                    dealloc_list,   mem)

# Language specific comment operations
comment(::LangMatlab,s)="% $s";
comment(::LangJulia,s)="# $s";
comment(::LangC,s)="/* $s */";

# Language specific memory operations
slotname(::LangJulia,i)="memslots[$i]"
slotname(::LangMatlab,i)="memslots{$i}"
slotname(::LangC,i)="memslots+($i-1)*n*n"

# Language specific variable declaration, initialization, and assignment
# For Julia more parsing may be required TODO
assign_coeff(::LangJulia,v,i)=("coeff$i","coeff$i=$v");
assign_coeff(::LangMatlab,v,i)=("coeff$i","coeff$i=$(real(v)) + 1i*$(imag(v))");
# In C, complex types are structures and are passed by reference.
function assign_coeff(::LangC_MKL,val::T,i) where T<:Complex
    assignment_string = "coeff$i.real = "*string(real(val))*";\n"*
                        "coeff$i.imag = "*string(imag(val))*";"
    variable_string = "&coeff$i"
    return (variable_string,assignment_string)
end
function assign_coeff(::LangC_OpenBLAS,val::T,i) where T<:Complex
    assignment_string = "coeff$i = "*string(real(val))*" + "*
                                     string(imag(val))*"*I;"
    variable_string = "&coeff$i"
    return (variable_string,assignment_string)
end
function assign_coeff(::LangC,val::T,i) where T<:Real
    assignment_string = "coeff$i = $val;";
    variable_string = "coeff$i"
    return (variable_string,assignment_string)
end

# Language specific constant declaration and reference
function declare_constant(::LangC_MKL,val::T,id,type) where T<:Complex
    return "const $type $id = {.real = "*string(real(val))*", "*
                              ".imag = "*string(imag(val))*"};"
end
function declare_constant(::LangC_OpenBLAS,val::T,id,type) where T<:Complex
    return "const $type $id = "*string(real(val))*" + "*
                                string(imag(val))*"*I;"
end
function declare_constant(::LangC,val::T,id,type) where T<:Real
    return "const $type $id = $val;"
end
reference_constant(::LangC,T,id)=(T<:Complex) ? "&$id" : "$id"





### Julia

function function_definition(lang::LangJulia,T,funname)
    code=init_code(lang);
    push_code!(code,"using LinearAlgebra",ind_lvl=0);
    push_code!(code,"function $funname(A)",ind_lvl=0);
    return code
end
function function_init(lang::LangJulia,T,mem,graph)
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

    push_code!(code,"memslots[j]=Matrix{T}(undef,n,n);",ind_lvl=2);
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
    push_code!(code,"end",ind_lvl=0);
    return code
end

# Adding of an identity matrix in julia code
function execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
    push_comment!(code,"Add lincomb with identity. Use view of diaganal.");
    push_code!(code,"copy!($(nodemem),$(non_I_parent_mem))");
    push_code!(code,"$(nodemem) .*= $non_I_parent_coeff");
    push_code!(code,"D=view($(nodemem), diagind($(nodemem), 0));");
    push_code!(code,"D .+= $I_parent_coeff");
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
    parent1mem=get_slot_name(mem,parent1)
    parent2mem=get_slot_name(mem,parent2)

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

            if (lang.exploit_uniformscaling &&  (parent1 == :I || parent2 == :I))
                # BLAS does not work with unform scaling use inplace instead
                if (parent1 == :I)
                    non_I_parent_mem=parent2mem;
                    non_I_parent_coeff=coeff2;
                    I_parent_coeff=coeff1;
                elseif (parent2 == :I)
                    non_I_parent_mem=parent1mem;
                    non_I_parent_coeff=coeff1;
                    I_parent_coeff=coeff2;
                end
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)



            else

                if (recycle_parent == parent1)
                    push_code!(code,"BLAS.axpby!($coeff2,$parent2mem,$coeff1,$nodemem);");
                else
                    push_code!(code,"BLAS.axpby!($coeff1,$parent1mem,$coeff2,$nodemem);");
                end
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

            if (parent1 == :I)

                non_I_parent_mem=parent2mem;
                non_I_parent_coeff=coeff2;
                I_parent_coeff=coeff1;
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
            elseif (parent2 == :I)
                non_I_parent_mem=parent1mem;
                non_I_parent_coeff=coeff1;
                I_parent_coeff=coeff2;
                execute_julia_I_op(code,nodemem,non_I_parent_mem,non_I_parent_coeff,I_parent_coeff)
            else
                push_code!(code,
                           "$(nodemem)[:]=$coeff1*$parent1mem + $coeff2*$parent2mem")

            end


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
function function_definition(lang::LangMatlab,T,funname)
    code=init_code(lang);
    push_code!(code,"function output=$funname(A)",ind_lvl=0);
    return code
end

function function_init(lang::LangMatlab,T,mem,graph)
    code=init_code(lang);
    push_code!(code,"n=size(A,1);");
    push_code!(code,"I=eye(n,n);");
    return code;
end

function init_mem(lang::LangMatlab,max_nof_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i));
    return mem;
end
function function_end(lang::LangMatlab,graph,mem)
    code=init_code(lang);
    push_code!(code,"output=$(graph.outputs[end]);");
    push_code!(code,"end",ind_lvl=0)
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
        push_code!(code,"$node= $coeff1*$parent1+$coeff2*$parent2;")
    end
    return (code,"$node");
end





### C

# Return MKL/OpenBLAS-specific includes.
function get_blas_includes(::LangC_MKL)
    return "#include<mkl/mkl.h>"
end
function get_blas_includes(::LangC_OpenBLAS)
    return "#include<cblas.h>\n#include<lapacke.h>"
end

# Return type and BLAS/LAPACK prefix corresponding to type T.
function get_blas_type(::LangC,T::Type{Float32})
    return ("float","s")
end
function get_blas_type(::LangC,T::Type{Float64})
    return ("double","d")
end
function get_blas_type(::LangC_MKL,T::Type{Complex{Float32}})
    return ("MKL_Complex8","c")
end
function get_blas_type(::LangC_MKL,T::Type{Complex{Float64}})
    return ("MKL_Complex16","z")
end
function get_blas_type(::LangC_OpenBLAS,T::Type{Complex{Float32}})
    return ("openblas_complex_float","c")
end
function get_blas_type(::LangC_OpenBLAS,T::Type{Complex{Float64}})
    return ("openblas_complex_double","z")
end

# Memory management functions.
function langc_get_slot_name(mem,key)
    return key==:A ? "A" : get_slot_name(mem,key)
end

function function_definition(lang::LangC,T,funname)
    (blas_type,blas_prefix)=get_blas_type(lang,T)
    code=init_code(lang);
    push_code!(code,get_blas_includes(lang),ind_lvl=0)
    push_code!(code,"#include<assert.h>",ind_lvl=0)
    push_code!(code,"#include<stdlib.h>",ind_lvl=0)
    push_code!(code,"#include<string.h>",ind_lvl=0)
    push_code!(code,"")
    push_comment!(code, "Code for polynomial evaluation.",ind_lvl=0)
    push_code!(code,"void $blas_prefix$funname(const $blas_type *A, "*
        "const size_t n, $blas_type *output) {",ind_lvl=0)
    return code
end

function init_mem(lang::LangC,max_nof_nodes;)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i));
    return mem;
end

function function_init(lang::LangC,T,mem,graph)
    (blas_type,blas_prefix)=get_blas_type(lang,T)
    code=init_code(lang)
    max_nodes=size(mem.slots,1)

    # Initialization.
    graph_ops=values(graph.operations)
    push_code!(code,"size_t max_memslots = $max_nodes;")
    push_code!(code,"")
    push_comment!(code,"Initializations.")
    # Allocate coefficients for linear combinations, if needed.
    if :lincomb in graph_ops
        push_code!(code,"$blas_type coeff1, coeff2;")
    end
    # Declare constants ZERO and ONE, if needed.
    if :mult in graph_ops
        push_code!(code, declare_constant(lang,convert(T,0.),"ZERO",blas_type))
    end
    if :lincomb in graph_ops || :mult in graph_ops
        push_code!(code, declare_constant(lang,convert(T,1.),"ONE",blas_type))
    end
    # Array ipiv for pivots of GEPP (needed only if graph has linear systems).
    if :ldiv in graph_ops
        push_code!(code, "lapack_int *ipiv = malloc(n*sizeof(*ipiv));")
    end
    if max_nodes > 0
        push_code!(code,"$blas_type *memslots = malloc(n*n*max_memslots"*
                        "*sizeof(*memslots));")
    end

    # Store identity explicitly only if graph has a linear combination of I.
    if has_identity_lincomb(graph)
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,:I);
        push_code!(code,"size_t j;")
        push_code!(code,"memset($nodemem, 0, n*n*sizeof(*memslots));")
        push_code!(code,"for(j=0; j<n*n; j+=n+1)");
        push_code!(code,"memslots[j] = ONE;",ind_lvl=2)
    end

    return code
end

function function_end(lang::LangC,graph,mem)
    code=init_code(lang);
    retval_node=graph.outputs[end];
    retval=langc_get_slot_name(mem,retval_node);
    push_code!(code,"");
    push_comment!(code,"Prepare output.")
    push_code!(code,"memcpy(output, $retval, n*n*sizeof(*output));")
    if :ldiv in values(graph.operations)
        push_code!(code,"free(ipiv);")
    end
    if size(mem.slots,1) > 0
        push_code!(code,"free(memslots);")
    end
    push_code!(code,"}",ind_lvl=0);
    return code
end

function execute_operation!(lang::LangC,T,graph,node,dealloc_list,mem)
    (blas_type,blas_prefix)=get_blas_type(lang,T)

    op=graph.operations[node]
    parent1=graph.parents[node][1]
    parent2=graph.parents[node][2]

    # Keep deallocation list (used for smart memory management).
    dealloc_list=deepcopy(dealloc_list)
    setdiff!(dealloc_list,[:A]) # A is not stored explicitly.
    setdiff!(dealloc_list,keys(mem.special_names))

    # Check if parents have a memory slot.
    if has_identity_lincomb(graph) # Linear combination of I.
        p1_is_identity=p2_is_identity=false
    else
        p1_is_identity=parent1==:I
        p2_is_identity=parent2==:I
        setdiff!(dealloc_list,[:I])
    end

    # Get memory slots of parents, if explicitly stored.
    if !p1_is_identity
        parent1mem=langc_get_slot_name(mem,parent1)
    end
    if !p2_is_identity
        parent2mem=langc_get_slot_name(mem,parent2)
    end

    rzero=reference_constant(lang,T,"ZERO")
    rone=reference_constant(lang,T,"ONE")

    code=init_code(lang);
    push_code!(code,"");
    push_comment!(code,"Computing $node with operation: $op");
    if op == :mult
        (nodemem_i,nodemem)=get_free_slot(mem)
        alloc_slot!(mem,nodemem_i,node);
        push_code!(code,"cblas_$blas_prefix"*"gemm(CblasColMajor, "*
                        "CblasNoTrans, CblasNoTrans, n, n, n,\n"*
                        "            $rone, $parent1mem, n, $parent2mem, n,\n"*
                        "            $rzero, $nodemem, n);")
    elseif op == :ldiv
        # Find a memory slot for the LU factors.
        # As ?detrs computes the LU decomposition of parent1 in-place, we need
        # to make a copy of $parent1mem, unless parent1 is on the deallocation
        # list.
        if (parent1 in dealloc_list)
            push_comment!(code,"Reusing memory of $parent1 for LU factors.")
            lhsmem=parent1mem
            lhsmem_i=get_slot_number(mem,parent1)
            # Remove parent 1 from the deallocation list, but only temporarily.
            # This is to make sure that the LU factors do not get overwritten by
            # the result of the sytem solve.
            setdiff!(dealloc_list,[parent1])
        else
            (lhsmem_i,lhsmem)=get_free_slot(mem)
            alloc_slot!(mem,lhsmem_i,:dummy)
            push_code!(code,"memcpy($lhsmem, $parent1mem, "*
                            "n*n*sizeof(*memslots));")
        end

        # Compute LU decomposition of parent1.
        push_code!(code,"LAPACKE_$blas_prefix"*"getrf(LAPACK_COL_MAJOR, n, n, "*
                        "$lhsmem, n, ipiv);")

        # Solve the linear systems. We have two cases:
        if parent2 == :I
            # 1. If parent2 is the identity, we use the function ?getri, which
            # overwrites the LU factors with the matrix inverse.

            # Compute matrix inverse.
            push_code!(code,"LAPACKE_$blas_prefix"*"getri(LAPACK_COL_MAJOR, "*
                        "n, $lhsmem, n, ipiv);")

            # The LU factors should not be deallocated in this case, but we have
            # to point the slot lhsmem to the current node.
            alloc_slot!(mem,lhsmem_i,node);
        else
            # 2. If parent2 is not the identity, we use the function ?getrs,
            # which # overwrites B in AX = B. In this case, we need to copy
            # $parent2mem to # $nodemem first, unless parent2 is in the
            # deallocation list.

            # Allocate slot for result.
            if (parent2 in dealloc_list)
                push_comment!(code,"Reusing memory of $parent2 for solution.")
                nodemem=langc_get_slot_name(mem,parent2)
                setdiff!(dealloc_list,[parent2])
                set_slot_number!(mem,get_slot_number(mem,parent2),node)
            else
                (nodemem_i,nodemem)=get_free_slot(mem)
                alloc_slot!(mem,nodemem_i,node)
                push_code!(code,"memcpy($nodemem, $parent2mem, "*
                "n*n*sizeof(*memslots));")
            end

            # Solve linear system.
            push_code!(code,"LAPACKE_$blas_prefix"*"getrs(LAPACK_COL_MAJOR, "*
                            "'N', n, n,\n"*
                            "               $lhsmem, n, ipiv,\n"*
                            "               $nodemem, n);")

            # Deallocate LU factors.
            # TODO It might make sense to keep them if another system with the
            # same left-hand side is still in the graph.
            free!(mem,lhsmem_i)
        end

    elseif op == :lincomb
        coeff_vars=Vector{String}(undef,2);

        # Coefficients of graph should have type T.
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1);
        push_code!(code,coeff1_code);
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2);
        push_code!(code,coeff2_code);

        # Reuse parent in deallocation list for output (saves memcopy).
        if ((parent1 in dealloc_list) || (parent2 in dealloc_list))

            recycle_parent=(parent1 in dealloc_list) ? parent1 : parent2
            push_comment!(code,"Smart lincomb recycle $recycle_parent");

            # Perform operation.
            nodemem=langc_get_slot_name(mem,recycle_parent)
            if (recycle_parent == parent1)
                if p2_is_identity
                    push_code!(code,"cblas_$blas_prefix"*"scal(n*n, $coeff1, "*
                                    "$nodemem, 1);")
                    push_code!(code,"cblas_$blas_prefix"*"axpby(n, "*
                                    "$coeff2, &ONE, 0,\n"*
                                    "             $rone, $nodemem, n+1);");
                else
                    push_code!(code,"cblas_$blas_prefix"*"axpby(n*n, "*
                                    "$coeff2, $parent2mem, 1,\n"*
                                    "             $coeff1, $nodemem, 1);");
                end
            else
                if p1_is_identity
                    push_code!(code,"cblas_$blas_prefix"*"scal(n*n, $coeff2, "*
                                    "$nodemem, 1);")
                    push_code!(code,"cblas_$blas_prefix"*"axpby(n, "*
                                    "$coeff1, &ONE, 0,\n"*
                                    "             $rone, $nodemem, n+1);");
                else
                    push_code!(code,"cblas_$blas_prefix"*"axpby(n*n, "*
                                    "$coeff1, $parent1mem, 1,\n"*
                                    "             $coeff2, $nodemem, 1);");
                end
            end

            # Remove output from deallocation list.
            setdiff!(dealloc_list,[recycle_parent]);

            # Update CodeMem.
            j=get_slot_number(mem,recycle_parent);
            set_slot_number!(mem,j,node);
        else
            # Allocate new slot for result.
            (nodemem_i,nodemem)=get_free_slot(mem)
            alloc_slot!(mem,nodemem_i,node);

            # Compute linear combination.
            if p1_is_identity
                push_code!(code,"memcpy($nodemem, $parent2mem, "*
                                "n*n*sizeof(*memslots));")
                push_code!(code,"cblas_$blas_prefix"*"scal(n*n, $coeff2, "*
                                "$nodemem, 1);")
                push_code!(code,"cblas_$blas_prefix"*"axpby(n, "*
                                "$coeff1, &ONE, 0,\n"*
                                "             $rone, $nodemem, n+1);");
            elseif p2_is_identity
                push_code!(code,"memcpy($nodemem, $parent1mem, "*
                                "n*n*sizeof(*memslots));")
                push_code!(code,"cblas_$blas_prefix"*"scal(n*n, $coeff1, "*
                                "$nodemem, 1);")
                push_code!(code,"cblas_$blas_prefix"*"axpby(n*n, "*
                                "$coeff2, &ONE, 0,\n"*
                                "             $rone, $nodemem, n+1);");
            else
                push_code!(code,"memcpy($nodemem, $parent2mem, "*
                                "n*n*sizeof(*memslots));")
                push_code!(code,"cblas_$blas_prefix"*"axpby(n*n, "*
                                "$coeff1, $parent1mem, 1,\n"*
                                "             $coeff2, $nodemem, 1);")
            end
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

    nodemem = 0
    return (code,nodemem)
end

function scalar_to_string(::LangC,z)
    return "$z"
end
function scalar_to_string(::LangC_OpenBLAS,z::T) where T<:Complex
    return "$(real(z)) + $(imag(z))*I"
end
function scalar_to_string(::LangC_MKL,z::T) where T<:Complex
    return "{$(real(z)), $(imag(z))}"
end
function matrix_to_string(lang::LangC,A)
    (m,n) = size(A)
    matrix_string = ""
    for i=1:m
        for j=1:n
            matrix_string *=
                scalar_to_string(lang,A[i,j]) *
                (j!=n ? ", " : "")
        end
        matrix_string *= (i!=n ? ",\n" : "")
    end
    return matrix_string
end

compilation_string(::LangC_OpenBLAS,fname)=
    "gcc -o main_compiled $fname -labials -llapacke"
compilation_string(::LangC_MKL,fname)=
    "gcc -o main_compiled $fname -lmkl_rt"

function gen_main(lang::LangC,T,fname,funname)
    (blas_type,blas_prefix)=get_blas_type(lang,T)
    code=init_code(lang);
    push_code!(code,"")
    push_code!(code,"")
    push_code!(code,"")
    push_code!(code,"typedef $blas_type blas_type;",ind_lvl=0)
    push_code!(code,"")

    # Main function starts here.
    push_comment!(code,"Code snippet that calls $blas_prefix$funname().",
                  ind_lvl=0)
    push_comment!(code,"With the GNU Compiler Collection, compile with:",
                  ind_lvl=0)
    push_comment!(code,compilation_string(lang,fname),ind_lvl=0)
    push_code!(code,"int main() {",ind_lvl=0)
    push_code!(code,"size_t i;")

    # Generate matrix.
    n=3 # Size of dummy matrix.
    push_code!(code,"size_t n = $n;")
    push_code!(code,"blas_type A[9] = {")
    A=randn(T,n,n)
    push_code!(code,matrix_to_string(lang,A))
    push_code!(code,"};")

    # Call polynomial evaluation function.
    push_code!(code,"$blas_prefix$funname(A,n,A);")
    push_code!(code,"}",ind_lvl=0)
    return code
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

Currently supported language: `LangJulia`, `LangMatlab`,
 `LangC_MKL`, `LangC_OpenBLAS`.

"""
function gen_code(fname,graph;
                  priohelp=Dict{Symbol,Float64}(),
                  lang=LangJulia(),
                  funname="dummy",
                  generate_main=false)

    if has_trivial_nodes(graph)
        error("Please run compress_graph!() on the graph first.")
    end
    T=eltype(eltype(typeof(graph.coeffs.vals)));

    if (fname isa String)
        file = open(abspath(fname), "w+")
    else
        # Lazy: Print out to stdout if no filename
        file=Base.stdout;
    end

    (order, can_be_deallocated, max_nof_slots) =
        get_topo_order(graph; priohelp=priohelp);

    # max_nof_slots is the path width which gives
    # a bound on the number memory slots needed.


    println(file,to_string(function_definition(lang,T,funname)));

    # We do a double sweep in order to determine exactly
    # how many memory slots are needed. The first sweep
    # we carry out all operations but store only the
    # maximum number of memory slots needed. The
    # second sweep generates the code.


    mem=init_mem(lang,max_nof_slots+3)
    function_init(lang,T,mem,graph)

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
    println(file,to_string(function_init(lang,T,mem,graph)));


    for (i,node) in enumerate(order)
        (exec_code,result_variable)=execute_operation!(lang,
                                                       T,graph,node,
                                                       can_be_deallocated[i],
                                                       mem)
        println(file,to_string(exec_code))
    end
    println(file,to_string(function_end(lang,graph,mem)));

    # Generate main function, if necessary.
    if generate_main && typeof(lang) <: LangC
        exec_code=gen_main(lang,T,fname,funname)
        println(file,to_string(exec_code))
    end

    if (fname isa String)
        close(file)
    end

end
