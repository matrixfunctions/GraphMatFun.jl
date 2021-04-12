export gen_exp_native_jl


"""
     (graph,crefs)=gen_exp_native_jl(A)

Creates a graph for the native scaling-and-squaring for the matrix exponential, as
implemented in Julia. The matrix `A` is taken as
input to determine the length of the Padé approximant, and the number of squares applied.

References:

* N. J. Higham. Functions of Matrices. SIAM publications, Philadelphia, PA, 2008.

* N. J. Higham. The Scaling and Squaring Method for the Matrix Exponential Revisited. SIAM J. Matrix Anal. Appl., 2005 26:4, 1179-1193

* Julia's matrix exponential, at the time of conversion: https://github.com/JuliaLang/julia/blob/697e782ab86bfcdd7fd15550241fe162c51d9f98/stdlib/LinearAlgebra/src/dense.jl#L554
    """
function gen_exp_native_jl(A)
    nA = opnorm(A, 1)
    T = eltype(A)
    graph = Compgraph(T)
    if (nA <= 2.1)
        return gen_exp_native_jl_A(graph, nA, T)
    else
        return gen_exp_native_jl_B(graph, nA, T)
    end
end


# For sufficiently small nA, use lower order Padé-Approximations
function gen_exp_native_jl_A(graph, nA, T)
    if nA > 0.95
        C = T[17643225600.,8821612800.,2075673600.,302702400.,
              30270240.,   2162160.,    110880.,     3960.,
              90.,         1.]
    elseif nA > 0.25
        C = T[17297280.,8648640.,1995840.,277200.,
              25200.,   1512.,     56.,     1.]
    elseif nA > 0.015
        C = T[30240.,15120.,3360.,
              420.,   30.,   1.]
    else
        C = T[120.,60.,12.,1.]
    end

    n = length(C)
    s = (div(n, 2) - 1)
    cref=Vector{Tuple{Symbol,Int}}(undef,n)

    add_mult!(graph,:A2,:A,:A)
    evenmon = Vector{Symbol}(undef,s+1)
    evenmon[1] = :I
    evenmon[2] = :A2
    for k = 2:s
        sym = Symbol("A$(2*k)")
        add_mult!(graph,sym,:A2,evenmon[k])
        evenmon[k+1] = sym
    end

    a = view(C,1:2:n-1)
    cref[1:s+1] = add_sum!(graph, :V, a, evenmon, :V)

    a = view(C,2:2:n)
    cref[s+2:n] = add_sum!(graph, :Ua, a, evenmon, :Ua)
    add_mult!(graph,:U,:Ua,:A)

    add_lincomb!(graph,:X,1.0,:V,1.0,:U)
    add_lincomb!(graph,:Z,1.0,:V,-1.0,:U)
    add_ldiv!(graph,:P,:Z,:X)

    add_output!(graph,:P)

    return (graph, cref)
end


# Full scaling and squaring
function gen_exp_native_jl_B(graph, nA, T)
    s  = log2(nA/5.4) # power of 2 later reversed by squaring
    C=:A
    if s > 0
        si = ceil(Int,s)
        γ = 1/convert(T,2^si)
        C=:C
        add_lincomb!(graph,C,γ,:A,0,:I)
    end
    CC = T[64764752532480000.,32382376266240000.,7771770303897600.,
           1187353796428800.,  129060195264000.,  10559470521600.,
           670442572800.,      33522128640.,      1323241920.,
           40840800.,           960960.,           16380.,
           182.,                1.]

    cref=Vector{Tuple{Symbol,Int}}(undef,14)

    add_mult!(graph,:A2,C,C)
    add_mult!(graph,:A4,:A2,:A2)
    add_mult!(graph,:A6,:A2,:A4)


    #U  = A * (A6 * (CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2) .+
    #          CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn)

    # Ub3= CC[14].*A6 .+ CC[12].*A4 .+ CC[10].*A2
    a = view(CC,10:2:14)
    cref[1:3] = add_sum!(graph, :Ub3, a, [:A2, :A4, :A6], :Ub)
    add_mult!(graph,:Ub,:Ub3,:A6)

    # Ua = CC[8].*A6 .+ CC[6].*A4 .+ CC[4].*A2 .+ CC[2].*Inn
    a = view(CC,2:2:8)
    cref[4:7] = add_sum!(graph, :Ua, a, [:I, :A2, :A4, :A6], :Ua)

    add_lincomb!(graph,:Uc,1.0,:Ub,1.0,:Ua)
    add_mult!(graph,:U,C,:Uc)


    # V  = A6 * (CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2) .+
    #            CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn

    # Vb3= CC[13].*A6 .+ CC[11].*A4 .+ CC[9].*A2
    a = view(CC,9:2:13)
    cref[8:10] = add_sum!(graph, :Vb3, a, [:A2, :A4, :A6], :Vb)
    add_mult!(graph,:Vb,:Vb3,:A6)

    # Va = CC[7].*A6 .+ CC[5].*A4 .+ CC[3].*A2 .+ CC[1].*Inn
    a = view(CC,1:2:7)
    cref[11:14] = add_sum!(graph, :Va, a, [:I, :A2, :A4, :A6], :Va)

    add_lincomb!(graph,:V,1.0,:Vb,1.0,:Va)

    add_lincomb!(graph,:X,1.0,:V, 1.0,:U)
    add_lincomb!(graph,:Z,1.0,:V,-1.0,:U)
    add_ldiv!(graph,:P,:Z,:X)

    Qtm1=:P
    Qt=:P
    if (s>0)
        for t=1:si
            Qt=Symbol("S"*string(t))
            if (t>1)
                Qtm1=Symbol("S"*string(t-1))
            end
            add_mult!(graph,Qt,Qtm1,Qtm1)
        end
    end

    add_output!(graph,Qt)

    return (graph,cref)
end
