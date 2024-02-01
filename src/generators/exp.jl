export graph_exp_native_jl, graph_exp_native_jl_degopt

function get_expm_pade_coeffs(d, T)
    if d == 13
        return T[
            64764752532480000.0,
            32382376266240000.0,
            7771770303897600.0,
            1187353796428800.0,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0,
        ]
    elseif d == 9
        C = T[
            17643225600.0,
            8821612800.0,
            2075673600.0,
            302702400.0,
            30270240.0,
            2162160.0,
            110880.0,
            3960.0,
            90.0,
            1.0,
        ]
    elseif d == 7
        C = T[
            17297280.0,
            8648640.0,
            1995840.0,
            277200.0,
            25200.0,
            1512.0,
            56.0,
            1.0,
        ]
    elseif d == 5
        C = T[30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0]
    elseif d == 3
        C = T[120.0, 60.0, 12.0, 1.0]
    else
        error("Unknown degree for expm-Padé coefficients.")
    end
end

"""
     (graph,crefs)=graph_exp_native_jl(A; input=:A)

Creates a graph for the native scaling-and-squaring for the matrix exponential,
as implemented in Julia. The matrix `A` is taken as input to determine the
length of the Padé approximant and the number of squares applied, as well as
determining the type of the coefficients. The kwarg `input` determines the name
of the matrix, in the graph.

**References**

1. N. J. Higham. "Functions of Matrices". SIAM, Philadelphia, PA, 2008.
   DOI: [10.1137/1.9780898717778](https://doi.org/10.1137/1.9780898717778)

2. N. J. Higham. "The Scaling and Squaring Method for the Matrix Exponential
   Revisited". SIAM Journal on Matrix Analysis and Applications,
   26(4):1179-1193, 2005.
   DOI: [10.1137/04061101X](https://doi.org/10.1137/04061101X)

3. "Julia's matrix exponential", [at the time of
   conversion](https://github.com/JuliaLang/julia/blob/697e782ab86bfcdd7fd15550241fe162c51d9f98/stdlib/LinearAlgebra/src/dense.jl#L554).
"""
function graph_exp_native_jl(A; input = :A)
    T = eltype(A)
    nA = opnorm(A, 1)
    d = get_expm_pade_degree_native_lj(nA)
    if d < 13
        C = get_expm_pade_coeffs(d, T)
        return graph_exp_native_jl_low(C, input)
    else
        C = get_expm_pade_coeffs(d, T)
        s = log2(nA / 5.4) # power of 2 later reversed by squaring
        return graph_exp_native_jl_high(C, s, input)
    end
end

function get_expm_pade_degree_native_lj(nA)
    if nA <= 2.1
        if nA > 0.95
            return 9
        elseif nA > 0.25
            return 7
        elseif nA > 0.015
            return 5
        else
            return 3
        end
    else
        return 13
    end
end

# For sufficiently small nA, use lower order Padé-Approximations
function graph_exp_native_jl_low(C, input)
    T = eltype(C)
    graph = Compgraph(T)

    n = length(C)
    s = (div(n, 2) - 1)
    cref = Vector{Tuple{Symbol,Int}}(undef, n)

    add_mult!(graph, :A2, input, input)
    evenmon = Vector{Symbol}(undef, s + 1)
    evenmon[1] = :I
    evenmon[2] = :A2
    for k = 2:s
        sym = Symbol("A$(2*k)")
        add_mult!(graph, sym, :A2, evenmon[k])
        evenmon[k+1] = sym
    end

    a = view(C, 1:2:n-1)
    cref[1:s+1] = add_lincomb!(graph, :V, a, evenmon)

    a = view(C, 2:2:n)
    cref[s+2:n] = add_lincomb!(graph, :Ua, a, evenmon)
    add_mult!(graph, :U, :Ua, input)

    add_lincomb!(graph, :X, 1.0, :V, 1.0, :U)
    add_lincomb!(graph, :Z, 1.0, :V, -1.0, :U)
    add_ldiv!(graph, :P, :Z, :X)

    add_output!(graph, :P)

    return (graph, cref)
end

# Full scaling and squaring
function graph_exp_native_jl_high(CC, s, input)
    T = eltype(CC)
    graph = Compgraph(T)

    C = input
    if s > 0
        si = ceil(Int, s)
        γ = 1 / convert(T, 2^si)
        C = :C
        add_lincomb!(graph, C, γ, input, 0, :I)
    end

    cref = Vector{Tuple{Symbol,Int}}(undef, 14)

    add_mult!(graph, :A2, C, C)
    add_mult!(graph, :A4, :A2, :A2)
    add_mult!(graph, :A6, :A2, :A4)

    #U  = A * (A6 * (CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2) .+
    #          CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn)

    # Ub3= CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2
    a = view(CC, 10:2:14)
    cref[1:3] = add_lincomb!(graph, :Ub3, a, [:A2, :A4, :A6])
    add_mult!(graph, :Ub, :Ub3, :A6)

    # Ua = CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn
    a = view(CC, 2:2:8)
    cref[4:7] = add_lincomb!(graph, :Ua, a, [:I, :A2, :A4, :A6])

    add_lincomb!(graph, :Uc, 1.0, :Ub, 1.0, :Ua)
    add_mult!(graph, :U, C, :Uc)

    # V  = A6 * (CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2) .+
    #            CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn

    # Vb3= CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2
    a = view(CC, 9:2:13)
    cref[8:10] = add_lincomb!(graph, :Vb3, a, [:A2, :A4, :A6])
    add_mult!(graph, :Vb, :Vb3, :A6)

    # Va = CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn
    a = view(CC, 1:2:7)
    cref[11:14] = add_lincomb!(graph, :Va, a, [:I, :A2, :A4, :A6])

    add_lincomb!(graph, :V, 1.0, :Vb, 1.0, :Va)

    add_lincomb!(graph, :X, 1.0, :V, 1.0, :U)
    add_lincomb!(graph, :Z, 1.0, :V, -1.0, :U)
    add_ldiv!(graph, :P, :Z, :X)

    Qtm1 = :P
    Qt = :P
    if (s > 0)
        for t = 1:si
            Qt = Symbol("S" * string(t))
            add_mult!(graph, Qt, Qtm1, Qtm1)

            Qtm1 = Qt
        end
    end

    add_output!(graph, Qt)

    return (graph, cref)
end

"""
     (graph,crefs)=graph_exp_native_jl_degopt(A; input=:A)

Same as [`graph_exp_native_jl`](@ref) but with calls to [`graph_degopt`](@ref)
for contruction of numerator and denominator polynomials.
"""
function graph_exp_native_jl_degopt(A; input = :A)
    T = eltype(A)
    nA = opnorm(A, 1)
    d = get_expm_pade_degree_native_lj(nA)
    if d < 13
        C = get_expm_pade_coeffs(d, T)
        return graph_exp_native_jl_low_degopt(C, input)
    else
        C = get_expm_pade_coeffs(d, T)
        s = log2(nA / 5.4) # power of 2 later reversed by squaring
        return graph_exp_native_jl_high_degopt(C, s, input)
    end
end

# For sufficiently small nA, use lower order Padé-Approximations
function graph_exp_native_jl_low_degopt(C, input)
    T = eltype(C)
    graph = Compgraph(T)

    n = length(C)
    s = (div(n, 2) - 1)

    xU = Vector{Tuple{Vector{T},Vector{T}}}(undef, s + 1)
    # Even monomials
    xU[1] = (vcat(zeros(T, 1), one(T)), vcat(zeros(T, 1), one(T)))
    for k = 2:s
        xU[k] =
            (vcat(zeros(T, 2), one(T), zeros(k - 2)), vcat(zeros(T, k), one(T)))
    end

    # U
    a = view(C, 4:2:n)
    xU[s+1] =
        (vcat(zeros(T, 1), one(T), zeros(T, s )), vcat(C[2], zero(T), a))
    zU = vcat(zeros(T, s + 2), one(T))
    (graphU, crefU) = graph_degopt(xU, zU, input = input)
    rename_node!(graphU, graphU.outputs[1], :U, crefU)

    # V
    xV = Vector{Tuple{Vector{T},Vector{T}}}(undef, s)
    xV[1:s] = xU[1:s]
    a = view(C, 3:2:n-1)
    zV = vcat(C[1], zero(T), a)
    (graphV, crefV) = graph_degopt(xV, zV, input = input)
    rename_node!(graphV, graphV.outputs[1], :V, crefV)

    graph = merge_graphs(
        graphU,
        graphV,
        prefix1 = "U",
        prefix2 = "V",
        skip_basic1 = true,
        skip_basic2 = true,
        cref1 = crefU,
        cref2 = crefV,
        input1 = input,
        input2 = input,
    )
    cref = vcat(crefU, crefV)
    empty!(graph.outputs)

    merge_rename_exp_native_jl_degopt!(graph, cref, s + 1)

    add_lincomb!(graph, :X, 1.0, :VV, 1.0, :UU)
    add_lincomb!(graph, :Z, 1.0, :VV, -1.0, :UU)
    add_ldiv!(graph, :P, :Z, :X)

    add_output!(graph, :P)

    return (graph, cref)
end

# Full scaling and squaring
function graph_exp_native_jl_high_degopt(CC, s, input)
    T = eltype(CC)
    graph = Compgraph(T)
    C = input
    if s > 0
        si = ceil(Int, s)
        C = :C
        # Actual scaling below, after graph is created
    end

    #U  = A * (A6 * (CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2) .+
    #          CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn)
    xU = Vector{Tuple{Vector{T},Vector{T}}}(undef, 5)

    # Even monomials
    xU[1] = (vcat(zeros(T, 1), one(T)), vcat(zeros(T, 1), one(T)))
    xU[2] = (vcat(zeros(T, 2), one(T)), vcat(zeros(T, 2), one(T)))
    xU[3] = (vcat(zeros(T, 2), one(T), zero(T)), vcat(zeros(T, 3), one(T)))

    # Up = A6 * (CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2)
    a = view(CC, 10:2:14)
    xU[4] = (vcat(zeros(T, 4), one(T)), vcat(zeros(T, 2), a))

    # U = A* (CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn + Up)
    a = view(CC, 4:2:8)
    xU[5] =
        (vcat(zero(T), one(T), zeros(T, 4)), vcat(CC[2], zero(T), a, one(T)))

    zU = vcat(zeros(T, 6), one(T))
    (graphU, crefU) = graph_degopt(xU, zU, input = C)
    rename_node!(graphU, graphU.outputs[1], :U, crefU)

    # V  = A6 * (CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2) .+
    #            CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn
    xV = Vector{Tuple{Vector{T},Vector{T}}}(undef, 4)

    # Even monomials
    xV[1:3] = xU[1:3]

    # Vp = A6 * ( CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2 )
    a = view(CC, 9:2:13)
    xV[4] = (vcat(zeros(T, 4), one(T)), vcat(zeros(T, 2), a))

    # V = Vp + CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn
    a = view(CC, 3:2:7)
    zV = vcat(CC[1], zero(T), a, one(T))
    (graphV, crefV) = graph_degopt(xV, zV, input = C)
    rename_node!(graphV, graphV.outputs[1], :V, crefV)

    graph = merge_graphs(
        graphU,
        graphV,
        prefix1 = "U",
        prefix2 = "V",
        skip_basic1 = true,
        skip_basic2 = true,
        cref1 = crefU,
        cref2 = crefV,
        input1 = C,
        input2 = C,
    )
    cref = vcat(crefU, crefV)
    empty!(graph.outputs)

    merge_rename_exp_native_jl_degopt!(graph, cref, 4)

    add_lincomb!(graph, :X, 1.0, :VV, 1.0, :UU)
    add_lincomb!(graph, :Z, 1.0, :VV, -1.0, :UU)
    add_ldiv!(graph, :P, :Z, :X)

    if (s > 0)
        γ = 1 / convert(T, 2^si)
        add_lincomb!(graph, C, γ, input, 0, :I)
        xS = Vector{Tuple{Vector{T},Vector{T}}}(undef, si)
        for i = 1:si
            xS[i] = (vcat(zeros(T, i), one(T)), vcat(zeros(T, i), one(T)))
        end
        zS = vcat(zeros(T, si + 1), one(T))
        (graphS, crefS) = graph_degopt(xS, zS, input = :P)

        graph = merge_graphs(
            graph,
            graphS,
            prefix1 = "",
            prefix2 = "S",
            skip_basic1 = true,
            skip_basic2 = true,
            cref1 = cref,
            cref2 = crefS,
            input1 = C,
            input2 = :P,
        )
        cref = vcat(cref, crefS)
    else
        add_output!(graph, :P)
    end

    return (graph, cref)
end

# Remove redundant creation of monomial basis, e.g., VB3 == UB3 and VBa4_3 == UBa4_3
function merge_rename_exp_native_jl_degopt!(graph, cref, M)
    for N in ["VB", "UB"]
        for k in reverse(2:M)
            rename_node!(graph, Symbol(string(N, "$k")), Symbol("B$(k)"), cref)
        end
    end
end
