export LangC_MKL, LangC_OpenBLAS

# Data structures for the language.
"""
     LangC_MKL(gen_main::Bool)

Code generation in C using the oneAPI MKL implementation of BLAS.
"""
struct LangC_MKL
    gen_main::Any
    overwrite_input::Any
end

"""
     LangC_OpenBLAS(gen_main::Bool)

Code generation in C using the OpenBLAS implementation of BLAS.
"""
struct LangC_OpenBLAS
    gen_main::Any
    overwrite_input::Any
end
LangC = Union{LangC_MKL,LangC_OpenBLAS}

# Generated constructors
# https://docs.julialang.org/en/v1/manual/metaprogramming/#Code-Generation
# to default to gen_main=false, overwrite_input=true
for L in [:LangC_MKL, :LangC_OpenBLAS]
    @eval $L() = $L(false, true)
    @eval $L(gen_main) = $L(gen_main, true)
end

# Language specific operations.
comment(::LangC, s) = "/* $s */"
slotname(::LangC, i) = "memslots[$(i-1)]"

# Variable declaration, initialization, and reference.
# In C, complex types are structures and are passed by reference.
function initialization_string(::LangC, val::T) where {T<:Real}
    return "$val"
end
function initialization_string(::LangC_MKL, val::T) where {T<:Complex}
    return "{.real = " * string(real(val)) * ", " *
           ".imag = " * string(imag(val)) * "}"
end
function initialization_string(::LangC_OpenBLAS, val::T) where {T<:Complex}
    return string(real(val)) * " + " * string(imag(val)) * "*I"
end
function declare_var(lang::LangC, val, id, type)
    return "$type $id = " * initialization_string(lang, val) * ";"
end
function declare_coeff(lang::LangC, val, id, type)
    variable_string = "coeff_$id"
    dec_init_string = declare_var(lang, val, variable_string, type)
    return (variable_string, dec_init_string)
end

# Constant declaration and reference.
function declare_const(lang::LangC, val, id, type)
    return "const " * declare_var(lang::LangC, val, id, type)
end
reference_value(::LangC, T, id) = (T <: Complex) ? "&$id" : "$id"

# Code generation.
# Return MKL/OpenBLAS-specific includes.
function get_blas_includes(::LangC_MKL)
    return "#include<mkl.h>"
end
function get_blas_includes(::LangC_OpenBLAS)
    return "#include<cblas.h>\n#include<lapacke.h>"
end

# Return type and BLAS/LAPACK prefix corresponding to type T.
function get_blas_type(::LangC, T::Type{Float32})
    return ("float", "s")
end
function get_blas_type(::LangC, T::Type{Float64})
    return ("double", "d")
end
function get_blas_type(::LangC_MKL, T::Type{Complex{Float32}})
    return ("MKL_Complex8", "c")
end
function get_blas_type(::LangC_MKL, T::Type{Complex{Float64}})
    return ("MKL_Complex16", "z")
end
function get_blas_type(::LangC_OpenBLAS, T::Type{Complex{Float32}})
    return ("openblas_complex_float", "c")
end
function get_blas_type(::LangC_OpenBLAS, T::Type{Complex{Float64}})
    return ("openblas_complex_double", "z")
end

# Add functions to work with MKL complex types.
function add_auxiliary_functions(code, ::LangC, T)
    return # In general, nothing is needed.
end

function add_auxiliary_functions(code, ::LangC_MKL, T::Type{Complex{Float32}})
    push_code!(code,
"void fma_MKL_Complex8(MKL_Complex8 *acc,
                       const MKL_Complex8 *a,
                       const MKL_Complex8 *b) {
    acc->real += a->real * b->real - a->imag * b->imag;
    acc->imag += a->real * b->imag + a->imag * b->real;
}

void acc_MKL_Complex8(MKL_Complex8 *acc,
                       const MKL_Complex16 *a) {
    acc->real += a->real;
    acc->imag += a->imag;
}", ind_lvl = 0)
end

function add_auxiliary_functions(code, ::LangC_MKL, T::Type{Complex{Float64}})
    push_code!(code,
"void fma_MKL_Complex16(MKL_Complex16 *acc,
                       const MKL_Complex16 *a,
                       const MKL_Complex16 *b) {
    acc->real += a->real * a->real - b->imag * b->imag;
    acc->imag += a->real * b->imag + a->imag * b->real;
}

void acc_MKL_Complex16(MKL_Complex16 *acc,
                       const MKL_Complex16 *a) {
    acc->real += a->real;
    acc->imag += a->imag;
}", ind_lvl = 0)
end

function preprocess_codegen(graph, lang::LangC)
    return graph
end

function function_definition(lang::LangC, graph, T, funname, precomputed_nodes)
    (blas_type, blas_prefix) = get_blas_type(lang, T)
    code = init_code(lang)
    push_code!(code, get_blas_includes(lang), ind_lvl = 0)
    push_code!(code, "#include<assert.h>", ind_lvl = 0)
    push_code!(code, "#include<stdlib.h>", ind_lvl = 0)
    push_code!(code, "#include<string.h>", ind_lvl = 0)
    push_code!(code, "")
    # Add auxiliary functions to work with MKL complex types in lincomb.
    if :lincomb ∈ values(graph.operations)
        add_auxiliary_functions(code, lang, T)
    end
    push_code!(code, "")
    push_comment!(code, "Code for matrix function evaluation.", ind_lvl = 0)
    input_variables =
        join(map(x -> "$blas_type *" * string(x), precomputed_nodes), ", ")
    push_code!(
        code,
        "void $blas_prefix$funname($input_variables, " *
        "const size_t n, $blas_type *output) {",
        ind_lvl = 0,
    )
    return code
end

function init_mem(lang::LangC, max_nof_nodes, precomputed_nodes)
    mem = CodeMem(max_nof_nodes, i -> slotname(lang, i))

    # Make sure the precomputed nodes have memslots
    for (i, n) in enumerate(precomputed_nodes)
        alloc_slot!(mem, i, n)
    end

    return mem
end

function function_init(lang::LangC, T, mem, graph, precomputed_nodes)
    (blas_t, blas_prefix) = get_blas_type(lang, T)
    code = init_code(lang)
    max_nodes = size(mem.slots, 1)

    # Initialization.
    graph_ops = values(graph.operations)
    # memalloc: matrices to allocate
    # memslots: numbers of pointers
    num_precomputed_nodes = size(precomputed_nodes, 1)
    num_memalloc =
        lang.overwrite_input ? max_nodes - num_precomputed_nodes : max_nodes
    num_preallocated_slots = max_nodes - num_memalloc
    push_code!(code, "size_t max_memalloc = $num_memalloc;")
    push_code!(code, "size_t max_memslots = $max_nodes;")
    push_code!(code, "")

    push_comment!(code, "Declarations and initializations.")
    # Declare constants ZERO and ONE, if needed.
    if :mult in graph_ops
        push_code!(code, declare_const(lang, convert(T, 0.0), "ZERO", blas_t))
        push_code!(code, declare_const(lang, convert(T, 1.0), "ONE", blas_t))
    end
    # Array ipiv for pivots of GEPP (needed only if graph has linear systems).
    if :ldiv in graph_ops
        push_code!(code, "lapack_int *ipiv = malloc(n * sizeof(*ipiv));")
    end
    push_code!(code, "size_t j;")
    push_code!(code, "")

    push_comment!(code, "Memory management.")
    push_code!(
        code,
        "$blas_t *master_mem = malloc(n * n * max_memalloc" *
        " * sizeof(*master_mem));",
    )
    push_code!(
        code,
        "$blas_t *memslots[max_memslots]; " *
        comment(lang, "As many slots as nodes"),
    )

    # Initialize the inputs
    for (i, n) in enumerate(precomputed_nodes)
        Ak_slot_name = get_slot_name(mem, n)
        if (lang.overwrite_input)
            # Just assign a pointer to the slot to allow overwrite
            push_code!(
                code,
                "$Ak_slot_name = $n; " * comment(lang, "Overwrite $n"),
            )
        else
            # Otherwise make a copy
            push_code!(
                code,
                "memcpy($Ak_slot_name, $n, n * n * sizeof(*master_mem));",
            )
        end
    end

    # Initialize pointers.
    push_comment!(code, "The other slots are pointers to allocated memory.")
    push_code!(
        code,
        "for (j = $(num_preallocated_slots-1); j < max_memalloc; j++)"
    )
    push_code!(code, "memslots[j+1] = master_mem + j * n * n;", ind_lvl = 2)

    # Store identity explicitly only if graph has a linear combination of I.
    if has_identity_lincomb(graph)
        push_comment!(code, "Graph has linear combination of identities.")
        push_comment!(code, "The matrix I is explicitly allocated.")
        alloc_slot!(mem, num_precomputed_nodes + 1, :I)
        nodemem = get_slot_name(mem, :I)
        push_code!(code, "memset($nodemem, 0, n * n * sizeof(*master_mem));")
        push_code!(code, "for(j = 0; j < n * n; j += n + 1)")
        push_code!(code, "*($nodemem+j) = ONE;", ind_lvl = 2)
    end

    return code
end

function function_end(lang::LangC, graph, mem)
    code = init_code(lang)
    retval_node = graph.outputs[end]
    retval = get_slot_name(mem, retval_node)
    push_code!(code, "")
    push_comment!(code, "Prepare output.")
    push_code!(code, "memcpy(output, $retval, n * n * sizeof(*output));")
    if :ldiv in values(graph.operations)
        push_code!(code, "free(ipiv);")
    end
    if size(mem.slots, 1) > 0
        push_code!(code, "free(master_mem);")
    end
    push_code!(code, "}", ind_lvl = 0)
    return code
end

function add_lincomb_body(code, lang::LangC, T, nodemem, coeff_names, parent_mems)
    println(T)
    rhs = join((coeff_names .* " * ") .* ("*(" .* parent_mems.* " + i)"), " + ")
    for_body = "*(" * nodemem * " + i) = " * (isempty(rhs) ? "0" : rhs) * ";"
    push_code!(code, for_body; ind_lvl = 2)
end

fma_function_name(T::Type{Complex{Float32}}) = "fma_MKL_Complex8"
fma_function_name(T::Type{Complex{Float64}}) = "fma_MKL_Complex16"
function add_lincomb_body(code, lang::LangC_MKL, T::Type{Complex{S}},
                          nodemem, coeff_names, parent_mems) where S <: Real
    if isempty(coeff_names)
        lhs = "*(" * nodemem * " + i)"
        push_code!(code, lhs * ".real = 0;", ind_lvl = 2)
        push_code!(code, lhs * ".imag = 0;", ind_lvl = 2)
    else
        for (it, (coeff, parent)) in enumerate(zip(coeff_names, parent_mems))
            statement = fma_function_name(T) *
                    "(" * nodemem * " + i, " *
                    reference_value(lang, T, coeff) * ", " *
                    parent * " + i);"
            push_code!(code, statement; ind_lvl = 2)
        end
    end
end

function add_lincomb_identity_body(code, lang::LangC, T, nodemem, coeff_id)
    push_code!(code, "*(" * nodemem * " + i) += " * coeff_id * ";", ind_lvl = 2)
end

acc_function_name(T::Type{Complex{Float32}}) = "acc_MKL_Complex8"
acc_function_name(T::Type{Complex{Float64}}) = "acc_MKL_Complex16"
function add_lincomb_identity_body(code, lang::LangC_MKL, T::Type{Complex{S}},
                                   nodemem, coeff_id) where S <: Real
    statement = acc_function_name(T) * "(" * nodemem * " + i, " *
        reference_value(lang, T, coeff_id) * ");"
    push_code!(code, statement, ind_lvl = 2)
end

function execute_operation!(lang::LangC, T, graph, node, dealloc_list, mem)
    (blas_type, blas_prefix) = get_blas_type(lang, T)

    op = graph.operations[node]

    code = init_code(lang)
    push_code!(code, "")
    push_comment!(code, "Computing $node with operation: $op")

    # Keep deallocation list (used for smart memory management).
    dealloc_list = deepcopy(dealloc_list)
    setdiff!(dealloc_list, keys(mem.special_names))

    rzero = reference_value(lang, T, "ZERO")
    rone = reference_value(lang, T, "ONE")

    if op == :mult
        # The graph has been compressed, thus we can assume that neither parent1
        # is the identity.
        parent1 = graph.parents[node][1]
        parent2 = graph.parents[node][2]
        parent1mem = get_slot_name(mem, parent1)
        parent2mem = get_slot_name(mem, parent2)
        (nodemem_i, nodemem) = get_free_slot(mem)
        alloc_slot!(mem, nodemem_i, node)
        push_code!(
            code,
            "cblas_$blas_prefix" *
            "gemm(CblasColMajor, " *
            "CblasNoTrans, CblasNoTrans, n, n, n,\n" *
            "            $rone, $parent1mem, n, $parent2mem, n,\n" *
            "            $rzero, $nodemem, n);",
        )
    elseif op == :ldiv
        # Left parent cannot be the identity.
        parent1 = graph.parents[node][1]
        parent1mem = get_slot_name(mem, parent1)
        # Parent 2 can be the identity, but this is dealt with below.
        parent2 = graph.parents[node][2]
        # Find a memory slot for the LU factors.
        # As ?detrs computes the LU decomposition of parent1 in-place, we need
        # to make a copy of $parent1mem, unless parent1 is on the deallocation
        # list.
        if (parent1 in dealloc_list)
            push_comment!(code, "Reusing memory of $parent1 for LU factors.")
            lhsmem = parent1mem
            lhsmem_i = get_slot_number(mem, parent1)
            # Remove parent 1 from the deallocation list, but only temporarily.
            # This is to make sure that the LU factors do not get overwritten by
            # the result of the sytem solve.
            setdiff!(dealloc_list, [parent1])
        else
            (lhsmem_i, lhsmem) = get_free_slot(mem)
            alloc_slot!(mem, lhsmem_i, :dummy)
            push_code!(
                code,
                "memcpy($lhsmem, $parent1mem, " * "n * n * sizeof(*master_mem));",
            )
        end

        # Compute LU decomposition of parent1.
        push_code!(
            code,
            "LAPACKE_$blas_prefix" *
            "getrf(LAPACK_COL_MAJOR, n, n, " *
            "$lhsmem, n, ipiv);",
        )

        # Solve the linear systems. We have two cases:
        if parent2 == :I
            # 1. If parent2 is the identity, we use the function ?getri, which
            # overwrites the LU factors with the matrix inverse.

            # Compute matrix inverse.
            push_code!(
                code,
                "LAPACKE_$blas_prefix" *
                "getri(LAPACK_COL_MAJOR, " *
                "n, $lhsmem, n, ipiv);",
            )

            # The LU factors should not be deallocated in this case, but we have
            # to point the slot lhsmem to the current node.
            alloc_slot!(mem, lhsmem_i, node)
        else
            # 2. If parent2 is not the identity, we use the function ?getrs,
            # which # overwrites B in AX = B. In this case, we need to copy
            # $parent2mem to # $nodemem first, unless parent2 is in the
            # deallocation list.

            # Allocate slot for result.
            if (parent2 in dealloc_list)
                push_comment!(code, "Reusing memory of $parent2 for solution.")
                nodemem = get_slot_name(mem, parent2)
                setdiff!(dealloc_list, [parent2])
                set_slot_number!(mem, get_slot_number(mem, parent2), node)
            else
                parent2mem = get_slot_name(mem, parent2)
                (nodemem_i, nodemem) = get_free_slot(mem)
                alloc_slot!(mem, nodemem_i, node)
                push_code!(
                    code,
                    "memcpy($nodemem, $parent2mem, " *
                    "n*n*sizeof(*master_mem));",
                )
            end

            # Solve linear system.
            push_code!(
                code,
                "LAPACKE_$blas_prefix" *
                "getrs(LAPACK_COL_MAJOR, " *
                "'N', n, n,\n" *
                "               $lhsmem, n, ipiv,\n" *
                "               $nodemem, n);",
            )

            # Deallocate LU factors.
            # TODO It might make sense to keep them if another system with the
            # same left-hand side is still in the graph.
            free!(mem, lhsmem_i)
        end

    elseif op == :lincomb

        setdiff!(dealloc_list, [:I])

        fused_sum = (join("x*" .* string.(graph.parents[node]), " + "))
        push_comment!(code, "Computing $node = $fused_sum")

        # Set coefficients.
        coeff_names = Vector()
        parent_mems = Vector()
        id_coefficient = 0
        for (i, v) in enumerate(graph.coeffs[node])
            n = graph.parents[node][i]
            if (n == :I)
                # Coefficient of identities.
                id_coefficient += v
            else
                # Coefficient of other nodes.
                (blas_t, blas_prefix) = get_blas_type(lang, T)
                (coeff_i, coeff_i_code) = declare_coeff(lang, v,
                                                 "$node" * "_" * "$i", blas_t)
                push_code!(code, coeff_i_code)
                push!(coeff_names, coeff_i)
                push!(parent_mems, get_slot_name(mem, n))
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
            push_comment!(code, "Smart lincomb recycle $recycle_parent.")
        end

        # Write the linear combination.
        push_code!(code, "for (size_t i = 0; i < n * n; i++) {")
        add_lincomb_body(code, lang, T, nodemem, coeff_names, parent_mems)
        push_code!(code, "}")
        push_code!(code, "")

        if id_coefficient != 0
            (blas_t, blas_prefix) = get_blas_type(lang, T)
            (coeff_id, coeff_id_code) = declare_coeff(lang,
                    id_coefficient, "$node" * "_0", blas_t)
            push_code!(code, coeff_id_code)
            # for statement in coeff_id_code
            #    push_code!(code, statement)
            # end
            push_code!(code, "for (size_t i = 0; i < n*n; i += n + 1) {")
            add_lincomb_identity_body(code, lang, T, nodemem, coeff_id)
            push_code!(code, "}")
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
        i = get_slot_number(mem, n)
        push_comment!(code, "Deallocating $n in slot $i")
        free!(mem, i)
    end

    nodemem = 0
    return (code, nodemem)
end

function scalar_to_string(::LangC, z)
    return "$z"
end
function scalar_to_string(::LangC_OpenBLAS, z::T) where {T<:Complex}
    return "$(real(z)) + $(imag(z)) * I"
end
function scalar_to_string(::LangC_MKL, z::T) where {T<:Complex}
    return "{$(real(z)), $(imag(z))}"
end
function print_indented_matrix(lang::LangC, code, A; ind_lvl = 1)
    (m, n) = size(A)
    for j = 1:n
        column_string = ""
        for i = 1:m
            column_string *=
                scalar_to_string(lang, A[i, j]) * (i != m ? ", " : "")
        end
        column_string *= (j != n ? "," : "")
        push_code!(code, column_string, ind_lvl = ind_lvl)
    end
end

function compilation_string(::LangC_OpenBLAS, fname)
    return "gcc -o main_compiled $fname -lblas -llapacke"
end
function compilation_string(::LangC_MKL, fname)
    return "gcc -o main_compiled $fname -lmkl_rt"
end

function gen_main(lang::LangC, T, fname, funname; A = 10::Union{Integer,Matrix})
    code = init_code(lang)
    if (lang.gen_main)
        (blas_type, blas_prefix) = get_blas_type(lang, T)
        push_code!(code, "")
        push_code!(code, "")
        push_code!(code, "")
        push_code!(code, "typedef $blas_type blas_type;", ind_lvl = 0)
        push_code!(code, "")

        # Main function starts here.
        push_comment!(
            code,
            "Code snippet that calls $blas_prefix$funname().",
            ind_lvl = 0,
        )
        push_comment!(
            code,
            "With the GNU Compiler Collection, compile with:",
            ind_lvl = 0,
        )
        push_comment!(code, compilation_string(lang, fname), ind_lvl = 0)
        push_code!(code, "int main() {", ind_lvl = 0)
        push_code!(code, "size_t i;")

        # Generate matrix.
        if isa(A, Matrix)
            n = LinearAlgebra.checksquare(A)
            push_code!(code, "size_t n = $n;")
            push_code!(code, "blas_type A[$(n*n)] = {")
            print_indented_matrix(lang, code, A, ind_lvl = 2)
            push_code!(code, "};")
        else # A is an integer
            n = A
            push_code!(code, "size_t n = $n;")
            push_code!(code, "srand(0);")
            push_code!(code, "blas_type *A = malloc(n * n * sizeof(*A));")
            push_code!(code, "for(i = 0; i < n * n; i++){")
            push_code!(code, "A[i] = rand() / (1.0 * RAND_MAX);", ind_lvl = 2)
            push_code!(code, "}")
        end

        # Call polynomial evaluation function.
        push_code!(code, "blas_type *B = malloc(n * n * sizeof(*A));")
        push_code!(code, "$blas_prefix$funname(A, n, B);")
        push_code!(code, "return 0;")
        push_code!(code, "}", ind_lvl = 0)
    end

    return code
end
