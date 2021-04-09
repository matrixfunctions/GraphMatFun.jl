using LinearAlgebra
@testset "Polynomial evaluation" begin

    A = [3 4 ; 5 6.6];

    coeff = collect(0.1:0.1:2.1)

    poly_gens = Dict("Monomial" => (:gen_monomial,:gen_monomial_recursive),
                     "Horner"   => (:gen_horner,:gen_horner_recursive),
                     "PS"       => (:gen_ps,:gen_ps_recursive) )

    for key = keys(poly_gens)
        @testset "$key" begin
            for n = 2:length(coeff)
                (graph_1,cref_1) = eval(poly_gens[key][1])(coeff[1:n])
                (graph_2,cref_2) = eval(poly_gens[key][2])(coeff[1:n])
                @test eval_graph(graph_2,A) ≈ eval_graph(graph_1,A)
                @test sum(values(graph_2.operations) .== :mult) == sum(values(graph_1.operations) .== :mult)
            end

            if key == "Horner" # Add test here if scaling is available
                s = -0.89
                (graph_1,cref_1) = eval(poly_gens[key][2])(coeff .* s.^(0:length(coeff)-1))
                (graph_2,cref_2) = eval(poly_gens[key][2])(coeff, scaling=s)
                @test eval_graph(graph_1,A) ≈ eval_graph(graph_2,A)
            end
        end
    end
end
