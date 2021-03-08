using LinearAlgebra
@testset "code gen" begin

    (graph,crefs)=gen_ps([3 4 0.9 0.11])
    A=[3 4.0; 5.5 0.1];
    lang=LangJulia(false,false);
    fname=tempname()*".jl";
    gen_code(fname,graph,lang=lang)
    include(fname);
    @test eval_graph(graph,A)≈dummy(A)
    rm(fname);
end
