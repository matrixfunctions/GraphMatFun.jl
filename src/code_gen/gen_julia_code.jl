export LangJulia

# Data structure for the language.
struct LangJulia
    overwrite_input::Any # Overwrite input
    inline::Any
    dot_fusing::Any # Allow dot fusion
end
"""
    LangJulia(overwrite_input=true,inline=true,dot_fusing=true)

Code generation in julia language, with optional overwriting of input, inlining
the function and optional usage of dot fusion.
"""
LangJulia() = LangJulia(true, true, true)
LangJulia(overwrite_input) = LangJulia(overwrite_input, true, true)
function LangJulia(overwrite_input, inline)
    return LangJulia(overwrite_input, inline, true)
end

# Language specific operations.
comment(::LangJulia, s) = "# $s"

slotname(::LangJulia, i) = "memslots$i"

function assign_coeff(lang::LangJulia, v, i)
    if v == 1
        return (
            "value_one",
            comment(lang, "Saving scalar multiplications using ValueOne()."),
        )
    else
        assign_coeff_basic(lang, v, i)
    end
end

function assign_coeff_basic(lang::LangJulia, v, i)
    T = typeof(v)
    if big(T) == T # High precision coeffs should be parsed as such
        if real(T) == T
            return ("coeff$i", "coeff$i=big\"$v\"")
        else
            vr = real(v)
            vi = imag(v)
            return ("coeff$i", "coeff$i=big\"$vr\" + big\"$vi\"*im")
        end
    else
        return ("coeff$i", "coeff$i=$v")
    end
end

function preprocess_codegen(graph, lang::LangJulia)
    if (lang.dot_fusing)
        return MultiLincombCompgraph(graph) # Merge many lincombs for dot fusion
    else
        return graph
    end
end

# Code generation.
function push_code_matfun_axpby_I!(code)
    # Add code for matfun_axpby! using UniformScaling.
    return push_code_verbatim_string!(
        code,
        "
# Compute X <- a X + b I.
function matfun_axpby!(X,a,b,Y::UniformScaling)
    m,n=size(X)
    if ~(a isa ValueOne)
        rmul!(X,a)
    end
    @inbounds for i=1:n
        X[i,i]+=(b isa ValueOne) ? 1 : b
    end
end\n",
    )
end

function push_code_matfun_axpby!(code)
    # Add code for generic matfun_axpby! function.
    return push_code_verbatim_string!(
        code,
        "
# Compute X <- a X + b Y.
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
end\n",
    )
end

function function_definition(
    lang::LangJulia,
    graph,
    T,
    funname,
    precomputed_nodes,
)
    code = init_code(lang)
    push_code!(code, "using LinearAlgebra", ind_lvl = 0)
    # If graph has linear combinations, add corresponding axpby functions.
    if any(values(graph.operations) .== :lincomb)
        push_code!(code, "\nstruct ValueOne; end\nValueOne()", ind_lvl = 0)
    end
    lincomb_nodes =
        filter(x -> graph.operations[x] == :lincomb, keys(graph.operations))
    lincomb_with_I =
        filter(y -> any(map(x -> x == :I, graph.parents[y])), lincomb_nodes)
    if (!isempty(lincomb_with_I))
        push_code_matfun_axpby_I!(code)
    end
    if (!isempty(setdiff!(lincomb_nodes, lincomb_with_I)))
        push_code_matfun_axpby!(code)
    end

    input_variables = join(precomputed_nodes, ", ")
    inline_string = lang.inline ? "@inline " : ""
    push_code!(
        code,
        inline_string * "function $funname($input_variables)",
        ind_lvl = 0,
    )
    # Generate version a bang version of the function.
    if (lang.overwrite_input)
        copy_input = "copy(" * join(precomputed_nodes, "), copy(") * ")"
        push_code!(code, "return $(funname)!($copy_input)")
        push_code!(code, "end", ind_lvl = 0)
        push_code!(code, "", ind_lvl = 0)
        push_code!(
            code,
            inline_string * "function $(funname)!($input_variables)",
            ind_lvl = 0,
        )
    end
    push_code!(
        code,
        "T=promote_type(eltype($(input_variables[1])),$T) " * comment(
            lang,
            "Make it work for many 'bigger' types (matrices and scalars)",
        ),
    )
    return code
end

function function_init(lang::LangJulia, T, mem, graph, precomputed_nodes)
    code = init_code(lang)
    max_nodes = size(mem.slots, 1)

    # Allocate memory for memory slots.
    push_code!(code, "max_memslots=$max_nodes")
    push_code!(code, "n=size(A,1)")
    # Allocate input.
    push_comment!(
        code,
        "The first slots are precomputed nodes $precomputed_nodes",
    )
    jj = size(precomputed_nodes,1);
    if (!lang.overwrite_input)
        jj = 0;
    end
    for i=jj+1:max_nodes
        thisslotname = slotname(lang, i);
        push_code!(code,"$thisslotname = similar(A,T)");
    end

    push_comment!(code,"Assign precomputed nodes memslots ");

    precomp_nodes_string = join(repeat("A",jj),",");
    for (i, n) in enumerate(precomputed_nodes)
        Ak_slot_name = get_slot_name(mem, n)
        if (lang.overwrite_input)
            # Just assign a pointer to the slot to allow overwrite
            push_code!(
                code,
                "$Ak_slot_name=$n " * comment(lang, "overwrite $n"),
            )
        else
            # Otherwise make a copy
            push_code!(code, "copy!($Ak_slot_name,$n)")
        end
    end


    # If needed, allocate identity matrix.
    if has_identity_lincomb(graph)
        push_comment!(code, "Graph has linear combination of identities.")
        push_comment!(code, "The matrix I is explicitly allocated.")
        alloc_slot!(mem, length(precomputed_nodes) + 1, :I)
        I_slot_name = get_slot_name(mem, :I)
        push_code!(code, "$I_slot_name=Matrix{Float64}(I,n,n)")
    else
        push_comment!(code, "Uniform scaling is exploited.")
        push_comment!(code, "No matrix I explicitly allocated.")
    end

    # If needed, initialize variable of type ValueOne for axpby with a=1 or b=1.
    if (any(map(x -> any(x .== 1), values(graph.coeffs))))
        push_code!(code, "value_one=ValueOne()")
    end

    return code
end

function init_mem(lang::LangJulia, max_nof_nodes, precomputed_nodes)
    mem = CodeMem(max_nof_nodes, i -> slotname(lang, i))

    # Make sure the precomputed nodes have memslots
    for (i, n) in enumerate(precomputed_nodes)
        alloc_slot!(mem, i, n)
    end

    return mem
end

function function_end(lang::LangJulia, graph, mem)
    code = init_code(lang)
    retval_node = graph.outputs[end]
    retval = get_slot_name(mem, retval_node)

    push_code!(
        code,
        "return $retval " * comment(lang, "Returning $retval_node"),
    )
    push_code!(code, "end", ind_lvl = 0)
    return code
end

# Add a scaled identity to a matrix in julia code.
function execute_julia_I_op(
    code,
    nodemem,
    non_I_parent_mem,
    non_I_parent_coeff,
    I_parent_coeff,
)
    push_comment!(code, "Add lincomb with identity.")

    if (nodemem != non_I_parent_mem)
        push_code!(code, "copy!($(nodemem),$(non_I_parent_mem))")
    else
        push_comment!(code, "No copy necessary. Inline identity multiple add")
    end

    return push_code!(
        code,
        "matfun_axpby!($(nodemem),$non_I_parent_coeff,$I_parent_coeff,I)",
    )
end

function execute_operation!(
    lang::LangJulia,
    T,
    graph::MultiLincombCompgraph,
    node,
    dealloc_list,
    mem,
)
    op = graph.operations[node]
    if (op != :lincomb)
        return execute_operation_basic!(lang, T, graph, node, dealloc_list, mem)
    else        # Multiple additions goes here

        # Don't overwrite the list
        dealloc_list = deepcopy(dealloc_list)
        setdiff!(dealloc_list, keys(mem.special_names))
        setdiff!(dealloc_list, [:I])

        code = init_code(lang)

        fused_sum = (join("x*" .* string.(graph.parents[node]), '+'))
        push_comment!(code, "Computing $node = $fused_sum")

        # Write the coeffs into appropriate vectors
        coeff_names = Vector()
        coeff_list = graph.coeffs[node]
        parent_mems = Vector()
        id_coeffs = Vector() # Vector of identity additions
        for (i, v) in enumerate(coeff_list)
            (coeff_i, coeff_i_code) = assign_coeff_basic(lang, v, i)
            push_code!(code, coeff_i_code)
            n = graph.parents[node][i]
            if (n != :I)
                push!(coeff_names, coeff_i)
                push!(parent_mems, get_slot_name(mem, n))
            else  # The identites are added in a different vector
                push!(id_coeffs, coeff_i)
            end
        end

        # Determine a target mem: nodemem
        nodemem = nothing
        recycle_parent = nothing
        for n in graph.parents[node]
            if (n in dealloc_list) # If it's about to deallocated
                nodemem = get_slot_name(mem, n)
                recycle_parent = n
                break
            end
        end
        if (isnothing(recycle_parent))
            (nodemem_i, nodemem) = get_free_slot(mem)
            alloc_slot!(mem, nodemem_i, node)
        else
            push_comment!(code, "Smart lincomb recycle $recycle_parent")
        end

        # Print the code
        rhs = join((coeff_names .* ".*") .* parent_mems, " .+ ")
        lhs = nodemem
        if !isempty(rhs)
            push_code!(code, lhs * " .= " * rhs)
        else
            push_code!(code, lhs * " .= 0")
        end

        # Adjust the result with inplace additions of identity
        for c in id_coeffs
            if VERSION > v"1.7.0-DEV.1240"
                push_code!(code, "mul!($nodemem,true,I*$c,true,true)")
            else
                push_comment!(
                    code,
                    "This julia version does not support inplace identity add.",
                )
                push_comment!(
                    code,
                    "Update to Julia 1.7.0 or newer for better performance.",
                )
                push_code!(code, "$nodemem+=I*$c")
            end
        end

        if !(isnothing(recycle_parent))
            # Avoid deallocate recycled parent
            setdiff!(dealloc_list, [recycle_parent])

            # Set the memory pointer
            j = get_slot_number(mem, recycle_parent)
            set_slot_number!(mem, j, node)
        end

        # Deallocate
        for n in dealloc_list
            if n != :I # No memory is ever allocated for the identity matrix.
                i = get_slot_number(mem, n)
                push_comment!(code, "Deallocating $n in slot $i")
                free!(mem, i)
            end
        end

        return (code, nodemem)
    end
end

function execute_operation!(lang::LangJulia, T, graph, node, dealloc_list, mem)
    return execute_operation_basic!(lang, T, graph, node, dealloc_list, mem)
end

# The general base case. Separated for dispatch.
function execute_operation_basic!(
    lang::LangJulia,
    T,
    graph,
    node,
    dealloc_list,
    mem,
)
    op = graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]

    # Don't overwrite the list
    dealloc_list = deepcopy(dealloc_list)
    setdiff!(dealloc_list, keys(mem.special_names))
    setdiff!(dealloc_list, [:I])
    if parent1 != :I || has_identity_lincomb(graph)
        parent1mem = get_slot_name(mem, parent1)
    end
    if parent2 != :I || has_identity_lincomb(graph)
        parent2mem = get_slot_name(mem, parent2)
    end

    code = init_code(lang)
    push_comment!(code, "Computing $node with operation: $op")
    if op == :mult
        # Multiplication has no inline
        (nodemem_i, nodemem) = get_free_slot(mem)
        alloc_slot!(mem, nodemem_i, node)

        push_code!(code, "mul!($nodemem,$parent1mem,$parent2mem)")

    elseif op == :ldiv
        # Left division
        (nodemem_i, nodemem) = get_free_slot(mem)
        alloc_slot!(mem, nodemem_i, node)
        if parent2 == :I
            push_code!(code, "$nodemem=inv($parent1mem)")
        else
            push_code!(code, "$nodemem .=$parent1mem\\$parent2mem")
        end

    elseif op == :lincomb
        coeff_vars = Vector{String}(undef, 2)
        (coeff1, coeff1_code) = assign_coeff(lang, graph.coeffs[node][1], 1)
        push_code!(code, coeff1_code)
        (coeff2, coeff2_code) = assign_coeff(lang, graph.coeffs[node][2], 2)
        push_code!(code, coeff2_code)

        # Use economical memory slots
        if ((parent1 in dealloc_list) || (parent2 in dealloc_list))
            # Smart / inplace: parentX can be used to store newly computed value

            # No allocation needed

            recycle_parent = (parent1 in dealloc_list) ? parent1 : parent2

            push_comment!(code, "Smart lincomb recycle $recycle_parent")

            nodemem = get_slot_name(mem, recycle_parent)

            if parent1 == :I || parent2 == :I
                # BLAS does not work with unform scaling use inplace instead
                if (parent1 == :I)
                    non_I_parent_mem = parent2mem
                    non_I_parent_coeff = coeff2
                    I_parent_coeff = coeff1
                elseif (parent2 == :I)
                    non_I_parent_mem = parent1mem
                    non_I_parent_coeff = coeff1
                    I_parent_coeff = coeff2
                end
                execute_julia_I_op(
                    code,
                    nodemem,
                    non_I_parent_mem,
                    non_I_parent_coeff,
                    I_parent_coeff,
                )
            else
                if (recycle_parent == parent1)
                    push_code!(
                        code,
                        "matfun_axpby!($nodemem,$coeff1,$coeff2,$parent2mem)",
                    )
                else
                    push_code!(
                        code,
                        "matfun_axpby!($nodemem,$coeff2,$coeff1,$parent1mem)",
                    )
                end
            end

            # Avoid deallocated
            setdiff!(dealloc_list, [recycle_parent])

            # Set the memory pointer
            j = get_slot_number(mem, recycle_parent)
            set_slot_number!(mem, j, node)
        else
            # Default behavior:  Allocate new slot
            (nodemem_i, nodemem) = get_free_slot(mem)

            alloc_slot!(mem, nodemem_i, node)

            if (parent1 == :I)
                non_I_parent_mem = parent2mem
                non_I_parent_coeff = coeff2
                I_parent_coeff = coeff1
                execute_julia_I_op(
                    code,
                    nodemem,
                    non_I_parent_mem,
                    non_I_parent_coeff,
                    I_parent_coeff,
                )
            elseif (parent2 == :I)
                non_I_parent_mem = parent1mem
                non_I_parent_coeff = coeff1
                I_parent_coeff = coeff2
                execute_julia_I_op(
                    code,
                    nodemem,
                    non_I_parent_mem,
                    non_I_parent_coeff,
                    I_parent_coeff,
                )
            else
                push_code!(code, "copy!($(nodemem),$parent1mem)")
                push_code!(
                    code,
                    "matfun_axpby!($(nodemem),$coeff1,$coeff2,$parent2mem)",
                )
            end
        end

    else
        error("Unknown operation")
    end

    # Deallocate
    for n in dealloc_list
        if n != :I # No memory is ever allocated for the identity matrix.
            i = get_slot_number(mem, n)
            push_comment!(code, "Deallocating $n in slot $i")
            free!(mem, i)
        end
    end

    return (code, nodemem)
end
