export LangJulia

# Data structure for the language.
struct LangJulia
    overwrite_input::Any # Overwrite input
    inline::Any
    alloc_function
    only_overwrite
end
default_alloc_function(k)="similar(A,T)"

"""
    LangJulia(overwrite_input=true,inline=true,alloc_function,only_overwrite=false)

Code generation in julia language, with optional overwriting of input. The
parameter `alloc_function` is a function of three parameters `alloc_function(k)`
where `k` is the memory slot (default is `alloc_function(k)=similar(A,T)`). The
`only_overwrite` specifies if `f` should be created if the overwrite funtion
    `f!` contains the actual code.
"""
LangJulia() = LangJulia(true, true, default_alloc_function, false)
LangJulia(overwrite_input) = LangJulia(overwrite_input, true, default_alloc_function, false)
function LangJulia(overwrite_input, inline)
    return LangJulia(overwrite_input, inline, default_alloc_function, false)
end

# Language specific operations.
comment(::LangJulia, s) = "# $s"

slotname(::LangJulia, i) = "memslots$i"

function assign_coeff(::LangJulia, v, i)
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

# Code generation.
function function_definition(
    lang::LangJulia,
    graph,
    T,
    funname,
    precomputed_nodes,
)
    code = init_code(lang)

    push_code!(code, "using LinearAlgebra", ind_lvl = 0)
    input_variables = join(precomputed_nodes, ", ")
    inline_string = lang.inline ? "@inline " : ""
    # Generate version a bang version of the function.
    if (lang.overwrite_input)
        if (!lang.only_overwrite)
            push_code!(
                code,
                inline_string * "function $funname($input_variables)",
                ind_lvl = 0,
            )
            #
            eltype_precomputed = join(map(x->"eltype($x)", precomputed_nodes),",");
            push_code!(code,"T=promote_type($eltype_precomputed,$T)");

            for n in precomputed_nodes
                push_code!(code,"$(n)_copy=similar($n,T); copyto!($(n)_copy, $n);");
            end
            copy_input = join(map(x->"$(x)_copy", precomputed_nodes),",");
            push_code!(code, "return $(funname)!($copy_input)")
            push_code!(code, "end", ind_lvl = 0)
        end

        push_code!(code, "", ind_lvl = 0)
        push_code!(
            code,
            inline_string * "function $(funname)!($input_variables)",
            ind_lvl = 0,
        )
    else
        push_code!(
            code,
            inline_string * "function $funname($input_variables)",
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
    push_comment!(code, "max_memslots=$max_nodes")
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
        push_code!(code,"$thisslotname = "*lang.alloc_function(i));
    end

    push_comment!(code,"Assign precomputed nodes memslots");

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

# The general base case. Separated for dispatch.
function execute_operation!(
    lang::LangJulia,
    T,
    graph,
    node,
    dealloc_list,
    mem,
)
    op = graph.operations[node]

    # Don't overwrite the list
    dealloc_list = deepcopy(dealloc_list)
    setdiff!(dealloc_list, keys(mem.special_names))
    setdiff!(dealloc_list, [:I])

    if op ∈ [:mult, :ldiv]
        parent1 = graph.parents[node][1]
        parent2 = graph.parents[node][2]
        if parent1 != :I
            parent1mem = get_slot_name(mem, parent1)
        end
        if parent2 != :I
            parent2mem = get_slot_name(mem, parent2)
        end
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
        fused_sum = (join("x*" .* string.(graph.parents[node]), '+'))
        push_comment!(code, "Computing $node = $fused_sum")

        # Write the coeffs into appropriate vectors
        coeff_names = Vector()
        coeff_list = graph.coeffs[node]
        parent_mems = Vector()
        id_coeffs = Vector() # Vector of identity additions
        for (i, v) in enumerate(coeff_list)
            (coeff_i, coeff_i_code) = assign_coeff(lang, v, i)
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
            push_code!(code, "mul!($nodemem, true, I*$c, true, true)")
        end

        if !(isnothing(recycle_parent))
            # Avoid deallocate recycled parent
            setdiff!(dealloc_list, [recycle_parent])

            # Set the memory pointer
            j = get_slot_number(mem, recycle_parent)
            set_slot_number!(mem, j, node)
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
