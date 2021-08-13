using LinearAlgebra, StaticArrays
@testset "code gen" begin
    a = Vector{Number}(undef, 6)
    a[1] = 0.11
    a[2] = 0.11 + 0.1im
    a[3] = Float32(0.11)
    a[4] = ComplexF32(0.11 + 0.1im)
    a[5] = big"0.11"
    a[6] = parse(Complex{BigFloat}, "0.11+0.1im")

    # Test julia code generation
    for i = 1:length(a)
        (graph, crefs) = graph_ps([3 4 2 a[i]])
        add_ldiv!(graph, :R0, :A2, :P0)
        clear_outputs!(graph)
        add_output!(graph, :R0)
        for t3 in (true, false)
            for t2 in (true, false)
                for t1 in (true, false)
                    @show i,t1,t2,t3
                    TT = eltype(a);
                    A = convert.(TT,[3 4.0; 5.5 0.1]);
                    lang = LangJulia(t1,t2,t3)
                    fname = tempname() * ".jl"
                    ti = time_ns()
                    gen_code(fname, graph, lang = lang, funname = "dummy_$(ti)")
                    # and execution
                    include(fname)
                    @test eval_graph(graph, A) ≈ eval(Symbol("dummy_$(ti)"))(A)
                    rm(fname)
                end
            end
        end

    end

    # Test Statically sized matrix
    (graph, crefs) = graph_ps([3 4 2 10.0])
    fname = tempname() * ".jl"
    gen_code(fname, graph, funname = "thisfunction")
    # and execution
    include(fname)
    rm(fname)

    A=[  0.704017  -0.358008   0.0649856;  -1.50413    0.0251634  0.0170085;  0.348764   0.334011   0.0915101]
    Am=MMatrix{3,3}(A);
    @test thisfunction(A)≈thisfunction(Am)

    # Test precomputed nodes
    fname = tempname() * ".jl"
    gen_code(fname, graph, funname = "thisfunction2",precomputed_nodes=[:A,:A2])
    # and execution
    include(fname)
    rm(fname)
    @test thisfunction2(Am,Am*Am)≈thisfunction(Am)

    for i = 1:4 #Not high-precision test for Matlab and C
        (graph, crefs) = graph_ps([3 4 2 a[i]])
        add_ldiv!(graph, :R0, :A2, :P0)
        clear_outputs!(graph)
        add_output!(graph, :R0)

        # Test Matlab code generation (not execution)
        fname = tempname() * ".m"
        gen_code(fname, graph, lang = LangMatlab())
        rm(fname)

        # Test C code generation
        fname = tempname() * ".c"
        gen_code(fname, graph, lang = LangC_MKL())
        rm(fname)

        fname = tempname() * ".c"
        gen_code(fname, graph, lang = LangC_MKL(true))
        rm(fname)

        fname = tempname() * ".c"
        gen_code(fname, graph, lang = LangC_OpenBLAS())
        rm(fname)

        fname = tempname() * ".c"
        gen_code(fname, graph, lang = LangC_OpenBLAS(true))
        rm(fname)
    end
end
