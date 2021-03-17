using LinearAlgebra
@testset "code gen" begin

    (graph,crefs)=gen_ps([3 4 0.9 0.11])
    A=[3 4.0; 5.5 0.1];

    # Test julia code generation
    lang=LangJulia(false,false);
    fname=tempname()*".jl";
    gen_code(fname,graph,lang=lang)
    # and execution
    include(fname);
    @test eval_graph(graph,A)≈dummy(A)
    rm(fname);

    # Test Matlab code generation (not execution)
    fname=tempname()*".m";
    gen_code(fname,graph,lang=LangMatlab())
    rm(fname);


end