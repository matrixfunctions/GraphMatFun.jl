using LinearAlgebra
@testset "Polynomial evaluation" begin

    A = [3 4 ; 5 6.6];

    coeff = [3.0; -1.0; 2.0; 0.1; 1.2];

    P1 = coeff[1]*I + coeff[2]*A + coeff[3]*A^2 + coeff[4]*A^3 + coeff[5]*A^4
    P2 = coeff[1]*I + coeff[2]*A
    P3 = coeff[1]*I

    poly_gens = Dict("Monomial" => :graph_monomial,
                     "Horner"   => :graph_horner,
                     "PS"       => :graph_ps)

    for key = keys(poly_gens)
        @testset "$key" begin
            (graph,cref) = eval(poly_gens[key])(coeff)
            @test eval_graph(graph,A) ≈ P1

            (graph,cref) = eval(poly_gens[key])(coeff[1:2])
            @test eval_graph(graph,A) ≈ P2

            (graph,cref) = eval(poly_gens[key])([coeff[1]])
            @test eval_graph(graph,A) ≈ P3

            if key == "Horner" # Add test here if scaling is available
                s = -0.89
                (graph_1,cref_1) = eval(poly_gens[key])(coeff .* s.^(0:length(coeff)-1))
                (graph_2,cref_2) = eval(poly_gens[key])(coeff, scaling=s)
                @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)
            end
        end
    end
end
