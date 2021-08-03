export graph_ps, graph_ps_degopt

"""
     (graph,crefs)=graph_ps(a; input=:A,
                          B_base=:B, C_base=:C, P_base=:P)

Generates the graph for the Paterson–Stockmeyer procedure with monomial basis
coefficieents. More precily, it corresponds to evaluation of the polynomial

    p(s) = a[1] + a[2]*s + ... + a[n]*s^(n-1).

The code follows the description in [^F19].

Reference:

[^PS73]: M. Paterson and L. Stockmeyer. "On the number of nonscalar
    multiplications necessary to evaluate polynomials". SIAM Journal on
    Scientific Computing, 2(1):60-66, 1973. DOI:
    [10.1137/0202007](https://doi.org/https://doi.org/)
[^F19]: M. Fasi. "Optimality of the Paterson–Stockmeyer method for evaluating
    matrix polynomials and rational matrix functions". Linear Algebra and its
    Applications, 574, 2019. DOI:
    [10.1016/j.laa.2019.04.001](https://doi.org/10.1016/j.laa.2019.04.001)
"""
function graph_ps(a; input = :A, B_base = :B, C_base = :C, P_base = :P)

    # Initial setup
    n = length(a)
    cref = Vector{Tuple{Symbol,Int}}(undef, n)
    graph = Compgraph(eltype(a))

    k = n - 1 # Polynomial degree

    # Degenerate cases k = 0 and k = 1
    if k < 2
        if k == 0
            add_lincomb!(graph, :P0, a[1], :I, 0.0, :I)
            cref[1] = (:P0, 1)
        elseif k == 1
            add_lincomb!(graph, :P0, a[2], :A, a[1], :I)
            cref[1] = (:P0, 1)
            cref[2] = (:P0, 2)
        end
        add_output!(graph, :P0)
        return (graph, cref)
    else
        s = ceil(Int, sqrt(k))
        v = floor(Int, k / s)
    end

    # Create monomial basis
    prevkey = input
    # up to s-1 are needed for sub-polynomials, and s needed for outer monomials
    for i = 2:s
        key = Symbol("$input$i")
        add_mult!(graph, key, input, prevkey)
        prevkey = key
    end
    Amax = Symbol("$input$s")

    # Create intermediate coefficient-polynomials. (3) and following equation in
    # Fasi 2019
    for i = 0:(v-1)
        # The number of coefficients in the sub-coefficient-polynomial
        # Note, degree is one lower than the index.
        lower_idx = i * s + 1 #plus 1 for 1-indexing
        upper_idx = lower_idx + s - 1 # s number of coeffients, so degree s-1
        idx = lower_idx:upper_idx
        c = view(a, idx)
        ccref = view(cref, idx)

        # B-keys are cumulative sums of the coefficient-polynomials.
        Bkey = Symbol("$(B_base)_$(i)_1")
        add_lincomb!(graph, Bkey, c[1], :I, c[2], input)
        ccref[1] = (Bkey, 1)
        ccref[2] = (Bkey, 2)
        for k = 3:s
            #As above, degree is one lower than the index
            powkey = Symbol("$input$(k-1)")
            prevBkey = Bkey
            Bkey = Symbol("$(B_base)_$(i)_$(k-1)")
            add_lincomb!(graph, Bkey, c[k], powkey, 1.0, prevBkey)
            ccref[k] = (Bkey, 1)
        end
    end

    # Handle highest powers separately. i=v. Part of (3), the following
    # equation, and (4) in Fasi 2019.
    ks = rem(k, s) # Degree of remaining highest order term. We compute B_v_ks
    if ks == 0 # Only constant term
        Pkey = Symbol("$(P_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v-1)_$(s-1)")
        # P_{v-1} = a[n]*A^s + B_[v-1}
        add_lincomb!(graph, Pkey, a[n], Amax, 1.0, Bkey)
        cref[n] = (Pkey, 1)
    else # ks > 0. At least constant and linear term.
        # Repeat above procedue with possibly lower-degree polynomial
        idx = (n-(ks+1)+1):n
        c = view(a, idx)
        ccref = view(cref, idx)

        # P_{v-1} = A^s*B_v + B_{v-1}, i.e., have to create B_v

        # B-keys are cumulative sums of the coefficient-polynomials.
        Bkey = Symbol("$(B_base)_$(v)_1")
        add_lincomb!(graph, Bkey, c[1], :I, c[2], input)
        ccref[1] = (Bkey, 1)
        ccref[2] = (Bkey, 2)
        for k = 3:ks+1
            # As above, degree is one lower than the index
            powkey = Symbol("$input$(k-1)")
            prevBkey = Bkey
            Bkey = Symbol("$(B_base)_$(v)_$(k-1)")
            add_lincomb!(graph, Bkey, c[k], powkey, 1.0, prevBkey)
            ccref[k] = (Bkey, 1)
        end
        # C_{v-1} = A^s * B_v
        Ckey = Symbol("$(C_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v)_$(ks)")
        add_mult!(graph, Ckey, Bkey, Amax)
        # P_{v-1} = C_{v-1}+ B_{v-1}
        Pkey = Symbol("$(P_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v-1)_$(s-1)")
        add_lincomb!(graph, Pkey, 1.0, Ckey, 1.0, Bkey)
    end

    prevPkey = Pkey

    # Evaluate outer-polynomial with Horner-type scheme. (4) in Fasi 2019.
    # Pout = :Nothing;
    Pout = Pkey
    for i in reverse(0:(v-2))
        # C_i = A^s * P_{i+1}
        Ckey = Symbol("$(C_base)$i")
        add_mult!(graph, Ckey, prevPkey, Amax)

        # P_i = C_i + B_i
        Bkey = Symbol("$(B_base)_$(i)_$(s-1)")
        Pkey = Symbol("$(P_base)$i")
        add_lincomb!(graph, Pkey, 1.0, Ckey, 1.0, Bkey)

        prevPkey = Pkey
        Pout = Pkey
    end

    add_output!(graph, Pout)
    return (graph, cref)
end

"""
     (graph,crefs)=graph_ps_degopt(a; input=:A)

Generates the same polynomial as [`graph_ps`](@ref), i.e., the graph for the
Paterson–Stockmeyer procedure with monomial basis coefficieents. However, it
does so by wrapping a call to [`graph_degopt`](@ref), resulting in more degrees
of freedom in `crefs`.
"""
function graph_ps_degopt(a; input = :A)

    # Initial setup
    n = length(a)
    T = eltype(a)
    graph = Compgraph(T)
    k = n - 1 # Polynomial degree

    if k <= 2 # Degenerate cases k = 0, k = 1, and k = 2
        if k == 0
            error("Does not implement degree-zero polynomial.")
        else
            # PS is not more efficient than monomial evaluation for low degrees,
            # but gives complications in code below.
            x = Vector{Tuple{Vector{T},Vector{T}}}(undef, k - 1)
            for i = 1:(k-1)
                x[i] = (
                    vcat(zero(T), one(T), zeros(i - 1)),
                    vcat(zeros(i), one(T)),
                )
            end
            return graph_degopt(x, a)
        end
    end

    s = ceil(Int, sqrt(k))
    v = floor(Int, k / s)
    ks = rem(k, s) # Degree of remaining highest order term.
    # Treat special if only constant. Save one multiplication
    if ks == 0
        i_adj = 1
    else
        i_adj = 0
    end

    x = Vector{Tuple{Vector{T},Vector{T}}}(undef, s + v - 1 - i_adj)
    # Create monomial basis (The first s-1 elements in the x-vector  <=>
    # The first s+1 slots in an x-component)
    for i = 1:(s-1)
        x[i] =
            (vcat(zeros(T, i), one(T)), vcat(zero(T), one(T), zeros(T, i - 1)))
    end

    # Evaluate outer-polynomial with Horner-type scheme, and at the same time
    # create the intermediate coefficient-polynomials
    if ks != 0 # If highest order term non-constant we create the polynomial
        idx = (n-(ks+1)+1):n
        diff = s - (ks + 1)
        c = view(a, idx)
        # A^s*P_{v-1}
        x[s] = (vcat(zeros(T, s), one(T)), vcat(c, zeros(T, diff), zero(T)))
    end

    for i in reverse(0:v-2)
        lower_idx = (i + 1) * s + 1
        upper_idx = lower_idx + s - 1
        idx = lower_idx:upper_idx
        c = view(a, idx)
        # Coeffs for A^s
        xl = vcat(zeros(T, s), one(T), zeros(T, v - 1 - i - i_adj))
        # If highest order term is constant we incorporate it here
        if (i == v - 2) && (ks == 0)
            xr = vcat(c, a[n]) # Coeffs for P_{v-i} + A^s*a[n]
        else
            # Coeffs for P_{v-i} + A^s*P_{v+1-i}
            xr = vcat(c, zero(T), zeros(T, v - 2 - i - i_adj), one(T))
        end
        x[(s-i_adj)+v-1-i] = (xl, xr) # A^s*( P_{v-i} + A^s*P_{v+1-i} )
    end

    z = vcat(a[1:s], zero(T), zeros(T, v - 1 - i_adj), one(T))
    return graph_degopt(x, z, input = input)
end
