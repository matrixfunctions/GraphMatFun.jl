using LinearAlgebra
@testset "monomial" begin

    A = [3 4 ; 5 6.6];

    coeff = [3; -1; 2.0; 0.1; 1.2];
    P = coeff[1]*I + coeff[2]*A + coeff[3]*A^2 + coeff[4]*A^3 + coeff[5]*A^4
    (graph,cref) = gen_monomial(coeff)
    @test eval_graph(graph,A) ≈ P

    coeff = [3; -1.0;];
    P = coeff[1]*I + coeff[2]*A
    (graph,cref) = gen_monomial(coeff)
    @test eval_graph(graph,A) ≈ P

    coeff = [3.0];
    P = coeff[1]*I
    (graph,cref) = gen_monomial(coeff)
    @test eval_graph(graph,A) ≈ P
end
