using LinearAlgebra, Polynomials
@testset "Polynomial evaluation" begin
    A = Polynomial("x")

    coeff = collect(1:21.0)

    poly_gens = Dict(
        "Monomial" => (:graph_monomial, :graph_monomial_degopt),
        "Horner"   => (:graph_horner, :graph_horner_degopt),
        "PS"       => (:graph_ps, :graph_ps_degopt),
    )

    for key in keys(poly_gens)
        @testset "$key" begin
            for n = 2:length(coeff)
                (graph_1, cref_1) = eval(poly_gens[key][1])(coeff[1:n])
                (graph_2, cref_2) = eval(poly_gens[key][2])(coeff[1:n])
                @test eval_graph(Compgraph(Any, graph_2), A) ==
                      eval_graph(Compgraph(Any, graph_1), A)
                @test sum(values(graph_2.operations) .== :mult) ==
                      sum(values(graph_1.operations) .== :mult)
            end

            if key == "Horner" # Add test here if scaling is available
                s = -0.89
                (graph_1, cref_1) =
                    eval(poly_gens[key][2])(coeff .* s .^ (0:length(coeff)-1))
                (graph_2, cref_2) = eval(poly_gens[key][2])(coeff, scaling = s)
                @test eval_graph(graph_1, A) ≈ eval_graph(graph_2, A)
            end
        end
    end
end
