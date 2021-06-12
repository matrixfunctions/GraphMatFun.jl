export LangMatlab

# Data structure for the language.
"""
    LangMatlab()

Code generation for the Matlab language.

"""
struct LangMatlab end

# Language specific operations
comment(::LangMatlab,s)="% $s"

slotname(::LangMatlab,i)="memslots{$i}"

function assign_coeff(::LangMatlab,v,i)
    if imag(v) == 0
        ("coeff$i","coeff$i=$(v)")
    else
        ("coeff$i","coeff$i=$(real(v)) + 1i*$(imag(v))")
    end
end



# Code generation.
function function_definition(lang::LangMatlab,graph,T,funname,precomputed_nodes)
    code=init_code(lang)
    input_variables=join(precomputed_nodes, ", ");
    push_code!(code,"function output=$funname($input_variables)",ind_lvl=0)
    return code
end

function function_init(lang::LangMatlab,T,mem,graph,precomputed_nodes)
    code=init_code(lang)
    push_code!(code,"n=size(A,1);")
    push_code!(code,"I=eye(n,n);")
    return code
end

function init_mem(lang::LangMatlab,max_nof_nodes,precomputed_nodes)
    mem=CodeMem(max_nof_nodes,i->slotname(lang,i))
    return mem
end

function function_end(lang::LangMatlab,graph,mem)
    code=init_code(lang)
    push_code!(code,"output=$(graph.outputs[end]);")
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

    code=init_code(lang)
    push_comment!(code,"Computing $node with operation: $op")
    # Needs to be updated
    if op == :mult
        push_code!(code,"$node=$parent1 * $parent2;")
    elseif op == :ldiv
        push_code!(code,"$node=$parent1 \\ $parent2;")
    elseif op == :lincomb
        (coeff1,coeff1_code)=assign_coeff(lang,graph.coeffs[node][1],1)
        push_code!(code,"$coeff1_code;")
        (coeff2,coeff2_code)=assign_coeff(lang,graph.coeffs[node][2],2)
        push_code!(code,"$coeff2_code;")
        push_code!(code,"$node= $coeff1*$parent1+$coeff2*$parent2;")
    end
    return (code,"$node")
end
