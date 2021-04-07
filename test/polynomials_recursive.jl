using LinearAlgebra
@testset "Polynomial evaluation" begin

    A = [3 4 ; 5 6.6];

    coeff = [3.0; -1.0; 2.0; 0.1; 1.2];

    poly_gens = Dict("Monomial" => (:gen_monomial,:gen_monomial_recursive),
                     "Horner"   => (:gen_horner,:gen_horner_recursive) )

    for key = keys(poly_gens)
        @testset "$key" begin
            (graph_1,cref_1) = eval(poly_gens[key][1])(coeff)
            (graph_2,cref_2) = eval(poly_gens[key][2])(coeff)
            @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)


            (graph_1,cref_1) = eval(poly_gens[key][1])(coeff[1:2])
            (graph_2,cref_2) = eval(poly_gens[key][2])(coeff[1:2])
            @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)

            if key == "Horner" # Add test here if scaling is available
                s = -0.89
                (graph_1,cref_1) = eval(poly_gens[key][2])(coeff .* s.^(0:length(coeff)-1))
                (graph_2,cref_2) = eval(poly_gens[key][2])(coeff, scaling=s)
                @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)
            end
        end
    end
end
