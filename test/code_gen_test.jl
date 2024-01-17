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
        (graph, crefs) = graph_ps_degopt([3 4 2 -1 2 a[i]])
        add_lincomb!(graph, :Q, [2.0], graph.outputs)# Check
        add_ldiv!(graph, :R2, :B4, :Q)
        add_lincomb!(graph, :Is, [2., 3., 5.], [:I, :I, :I])
        add_lincomb!(graph, :R1, [2., 2.], [:Is, :R2])
        add_ldiv!(graph, :R0, :R1, :I)
        clear_outputs!(graph)
        add_output!(graph, :R0)

        for t2 in (true, false)
            for t1 in (true, false)
                TT = eltype(a);
                A = convert.(TT,[3 4.0; 5.5 0.1]);
                ti = time_ns()
                lang = LangJulia(t1, t2)
                fname = tempname() * ".jl"
                begin # To avoid generated codes interfere
                    gen_code(fname, graph, lang = lang, funname = "dummy_$(ti)")
                    # and execution
                    include(fname)
                    @test eval_graph(graph, A) ≈ eval(Symbol("dummy_$(ti)"))(A)
                end

                rm(fname)
            end
        end

    end

    # Test Statically sized matrix
    (graph, crefs) = graph_ps_degopt([3 4 2 1 0.1 1.0])
    fname = tempname() * ".jl"
    lang = LangJulia(true, true)
    begin
        gen_code(fname, graph, funname = "thisfunction",lang=lang)
        # and execution
        include(fname)
        rm(fname)

        A=[  0.704017  -0.358008   0.0649856;  -1.50413    0.0251634  0.0170085;  0.348764   0.334011   0.0915101]
        Am=MMatrix{3,3}(A);
        @test thisfunction(A)≈thisfunction(Am)
    end

    # Test precomputed nodes
    fname = tempname() * ".jl"
    lang = LangJulia(true, true)
    begin
        gen_code(fname, graph, funname = "thisfunction2",precomputed_nodes=[:A,:A2])
        # and execution
        include(fname)
        rm(fname)
        @test thisfunction2(Am,Am*Am)≈thisfunction(Am)
    end

    for i = 1:4 #Not high-precision test for Matlab and C
        (graph, crefs) = graph_ps_degopt([3 4 2 a[i] 1 0])
        add_ldiv!(graph, :R1, :B4, graph.outputs[1])
        add_ldiv!(graph, :R2, :R1, :I)
        add_ldiv!(graph, :R3, :B4, :A)
        add_lincomb!(graph,:R4, [1., 2., 3.], [:I, :I, :I])
        add_lincomb!(graph,:R5, [2., 3., 4., 5., 6.], [:I, :R1, :R2, :R2, :R4])
        clear_outputs!(graph)
        add_output!(graph, :R5)

        # Test Matlab code generation (not execution)
        fname = tempname() * ".m"
        gen_code(fname, graph, lang = LangMatlab())
        rm(fname)

        # Test C code generation
        for gen_main in (true, false)
            for overwrite_input in (true, false)
                fname = tempname() * ".c"
                gen_code(fname, graph, lang = LangC_MKL(gen_main, overwrite_input))
                rm(fname)

                fname = tempname() * ".c"
                gen_code(fname, graph, lang = LangC_OpenBLAS(gen_main, overwrite_input))
                rm(fname)
            end
        end
    end
end
