using LinearAlgebra
@testset "code gen" begin

    a = Vector{Number}(undef,4)
    a[1] = 0.11; a[2] = 0.11+0.1im; a[3] = big"0.11"; a[4] = parse(Complex{BigFloat},"0.11+0.1im")
    A=[3 4.0; 5.5 0.1];

    # Test julia code generation
    for i = 1:length(a)
        (graph,crefs)=graph_ps([3 4 2 a[i]])
        for t1 in (true,false)
            lang=LangJulia(t1);
            fname=tempname()*".jl";
            gen_code(fname,graph,lang=lang)
            # and execution
            include(fname);
            @test eval_graph(graph,A)≈dummy(A)
            rm(fname);
        end
    end

    (graph,crefs)=graph_ps([3 4 2 a[1]])
    # Test Matlab code generation (not execution)
    fname=tempname()*".m";
    gen_code(fname,graph,lang=LangMatlab())
    rm(fname);

    # Test C code generation

    fname=tempname()*".c";
    gen_code(fname,graph,lang=LangC_MKL())
    rm(fname);

    fname=tempname()*".c";
    gen_code(fname,graph,lang=LangC_MKL(true))
    rm(fname);

    fname=tempname()*".c";
    gen_code(fname,graph,lang=LangC_OpenBLAS())
    rm(fname);

    fname=tempname()*".c";
    gen_code(fname,graph,lang=LangC_OpenBLAS(true))
    rm(fname);





end
