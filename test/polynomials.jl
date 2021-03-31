using LinearAlgebra
@testset "Polynomial evaluation" begin

    A = [3 4 ; 5 6.6];

    coeff = [3.0; -1.0; 2.0; 0.1; 1.2];
    s = -0.89;
    B = s*A;

    P1 = coeff[1]*I + coeff[2]*B + coeff[3]*B^2 + coeff[4]*B^3 + coeff[5]*B^4
    P2 = coeff[1]*I + coeff[2]*B
    P3 = coeff[1]*I

    poly_gens = Dict("Monomial" => :gen_monomial,
                     "Horner"   => :gen_horner,
                     "PS"       => :gen_ps)

    for key = keys(poly_gens)
        @testset "$key" begin
            (graph,cref) = eval(poly_gens[key])(coeff, scaling=s)
            @test eval_graph(graph,A) ≈ P1

            (graph,cref) = eval(poly_gens[key])(coeff[1:2], scaling=s)
            @test eval_graph(graph,A) ≈ P2

            (graph,cref) = eval(poly_gens[key])([coeff[1]], scaling=s)
            @test eval_graph(graph,A) ≈ P3
        end
    end
end
