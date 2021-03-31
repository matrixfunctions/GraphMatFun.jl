export gen_ps

"""
     (graph,crefs)=gen_ps(a; input=:A, scaling=1.0,
                          B_base=:B, C_base=:C, P_base=:P)

Generates the graph for the Paterson-Stockmayer procedure with monomial basis coefficieents. More precily, it corresponds to evaluation of the polynomial

    p(s) = a[1] + a[2]*(αs) + ... + a[n]*(αs)^(n-1)

where α=`scaling`.

Reference: The code is follows the description in

M. Fasi, Optimality of the Paterson–Stockmeyer method for evaluating matrix polynomials and rational matrix functions, Linear Algebra Appl. (2019)

    """
function gen_ps(a; input=:A, scaling=1.0,
                B_base=:B, C_base=:C, P_base=:P);

    # Initial setup
    n = length(a);
    cref=Vector{Tuple{Symbol,Int}}(undef,n)
    graph=Compgraph(eltype(a));

    k = n-1; # Polynomial degree

    if scaling != 1
        a = a .* scaling.^(0:k)
    end

    # Degenerate cases k = 0 and k = 1
    if k < 2
        if k == 0
            add_lincomb!(graph,:P0,a[1],:I,0.0,:I)
            cref[1] = (:P0,1);
        elseif k == 1
            add_lincomb!(graph,:P0,a[2],:A,a[1],:I)
            cref[1] = (:P0,1);
            cref[2] = (:P0,2);
        end
        add_output!(graph, :P0)
        return (graph,cref)
    else
        s = ceil(Int,sqrt(k));
        v = floor(Int,k/s);
    end

    # Create monomial basis
    prevkey = input;
    for i = 2:s # up to s-1 needed for sub-polynomials, and s needed for outer monomials
        key  = Symbol("$input$i");
        add_mult!(graph,key, input, prevkey);
        prevkey = key;
    end
    Amax = Symbol("$input$s")


    # Create intermediate coefficient-polynomials. (3) and following equation in Fasi 2019
    for i = 0:(v-1)
        # The number of coefficients in the sub-coefficient-polynomial
        # Note, degree is one lower than the index.
        lower_idx = i*s+1; #plus 1 for 1-indexing
        upper_idx = lower_idx+s-1; # s number of coeffients, so degree s-1
        idx = lower_idx:upper_idx
        c = view(a, idx)
        ccref = view(cref,idx)

        Bkey = Symbol("$(B_base)_$(i)_1") # B-keys are cumulative sums of the coefficient-polynomials.
        add_lincomb!(graph, Bkey, c[1], :I, c[2], input)
        ccref[1] = (Bkey,1);
        ccref[2] = (Bkey,2);
        for k = 3:s
            powkey = Symbol("$input$(k-1)"); #As above, degree is one lower than the index
            prevBkey = Bkey;
            Bkey = Symbol("$(B_base)_$(i)_$(k-1)");
            add_lincomb!(graph, Bkey, c[k], powkey, 1.0, prevBkey)
            ccref[k] = (Bkey,1);
        end
    end

    # Handle highest powers separately. i=v. Part of (3), the following equation, and (4) in Fasi 2019.
    ks = rem(k,s); # Degree of remaining highest order term. We compute B_v_ks
    if ks == 0 # Only constant term
        Pkey = Symbol("$(P_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v-1)_$(s-1)");
        # P_{v-1} = a[n]*A^s + B_[v-1}
        add_lincomb!(graph, Pkey, a[n], Amax, 1.0, Bkey)
        cref[n] = (Pkey,1)
    else # ks > 0. At least constant and linear term. Repeat above procedue with possibly lower-degree polynomial
        idx = (n-(ks+1)+1):n;
        c = view(a, idx)
        ccref = view(cref,idx)

        # P_{v-1} = A^s*B_v + B_{v-1}, i.e., have to create B_v
        Bkey = Symbol("$(B_base)_$(v)_1") # B-keys are cumulative sums of the coefficient-polynomials.
        add_lincomb!(graph, Bkey, c[1], :I, c[2], input)
        ccref[1] = (Bkey,1);
        ccref[2] = (Bkey,2);
        for k = 3:ks+1
            powkey = Symbol("$input$(k-1)"); #As above, degree is one lower than the index
            prevBkey = Bkey;
            Bkey = Symbol("$(B_base)_$(v)_$(k-1)");
            add_lincomb!(graph, Bkey, c[k], powkey, 1.0, prevBkey)
            ccref[k] = (Bkey,1);
        end
        # C_{v-1} = A^s * B_v
        Ckey = Symbol("$(C_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v)_$(ks)");
        add_mult!(graph, Ckey, Bkey, Amax)
        # P_{v-1} = C_{v-1}+ B_{v-1}
        Pkey = Symbol("$(P_base)$(v-1)")
        Bkey = Symbol("$(B_base)_$(v-1)_$(s-1)");
        add_lincomb!(graph, Pkey, 1.0, Ckey, 1.0, Bkey)
    end

    prevPkey = Pkey;

    # Evaluate outer-polynomial with Horner-type scheme. (4) in Fasi 2019.
    # Pout = :Nothing;
    Pout = Pkey;
    for i = reverse(0:(v-2))
        # C_i = A^s * P_{i+1}
        Ckey = Symbol("$(C_base)$i")
        add_mult!(graph, Ckey, prevPkey, Amax)

        # P_i = C_i + B_i
        Bkey = Symbol("$(B_base)_$(i)_$(s-1)");
        Pkey = Symbol("$(P_base)$i")
        add_lincomb!(graph, Pkey, 1.0, Ckey, 1.0, Bkey)

        prevPkey = Pkey;
        Pout = Pkey;
    end

    add_output!(graph, Pout)
    return (graph,cref)

end
