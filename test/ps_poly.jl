using LinearAlgebra
@testset "ps" begin

    A=[3 4 ; 5 6.6];
    x=[3; -1; 2.0; 0.1];
    P=x[1]*I+x[2]*A+x[3]*A^2+x[4]*A^3
    (graph,cref)=gen_ps(x)

    @test eval_graph(graph,A)≈ P

    coeff = [3; -1.0;];
    P = coeff[1]*I + coeff[2]*A
    (graph,cref) = gen_monomial(coeff)
    @test eval_graph(graph,A) ≈ P

    coeff = [3.0];
    P = coeff[1]*I
    (graph,cref) = gen_monomial(coeff)
    @test eval_graph(graph,A) ≈ P
end
