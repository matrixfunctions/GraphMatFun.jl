using LinearAlgebra
@testset "Polynomial evaluation" begin

    A = [3 4 ; 5 6.6];

    coeff = [3.0; -1.0; 2.0; 0.1; 1.2];
    s = -0.89;


    poly_gens = Dict("Monomial" => (:gen_monomial,:gen_monomial_recursive))

    for key = keys(poly_gens)
        @testset "$key" begin
            (graph_1,cref_1) = eval(poly_gens[key][1])(coeff, scaling=s)
            (graph_2,cref_2) = eval(poly_gens[key][2])(coeff, scaling=s)
            @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)


            (graph_1,cref_1) = eval(poly_gens[key][1])(coeff[1:2], scaling=s)
            (graph_2,cref_2) = eval(poly_gens[key][2])(coeff[1:2], scaling=s)
            @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)

        end
    end
end
