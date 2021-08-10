using LinearAlgebra
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
                    A = [3 4.0; 5.5 0.1]
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
