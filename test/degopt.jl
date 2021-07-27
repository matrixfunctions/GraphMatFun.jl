using LinearAlgebra
@testset "degopt" begin
    x1a = [3; 4.0]
    x1b = [0.1; 0.3]
    x2a = [0.1; 0.1; -0.1]
    x2b = [0.3; 0; -0.1]
    x = [(x1a, x1b), (x2a, x2b)]

    z = [1.0; -0.1; 0.1; 0.01]

    degopt1 = Degopt(deepcopy(x), deepcopy(z))
    # testing graph_degopt
    (g1, _) = graph_degopt(degopt1)

    # Testing grow
    degopt2 = grow!(deepcopy(degopt1))
    (g2, _) = graph_degopt(degopt1)
    @test eval_graph(g1, 0.3) == eval_graph(g2, 0.3)

    # Test get_degopt_crefs
    (x_crefs, y_crefs) = get_degopt_crefs(2)
    set_coeffs!(g1, 5.0, [x_crefs[1][2][2]])

    x[1][2][2] = 5.0
    degopt3 = Degopt(deepcopy(x), deepcopy(z))
    (g3, _) = graph_degopt(degopt3)

    @test eval_graph(g1, -3.0) == eval_graph(g3, -3.0)

    # Test scale!

    degopt4 = deepcopy(degopt1)
    scale!(degopt4, 2)
    (g4, _) = graph_degopt(degopt4)
    (g1, _) = graph_degopt(degopt1)
    @test eval_graph(g4, 0.3) == eval_graph(g1, 2 * 0.3)

    # Test square!
    degopt5 = deepcopy(degopt1)
    square!(degopt5)
    (g5, _) = graph_degopt(degopt5)
    @test eval_graph(g5, 2 + 2im) == eval_graph(g1, 2 + 2im)^2

    # Test initialize as degopt
    degopt6 = Degopt(g1)
    (g6, _) = graph_degopt(degopt6)
    @test eval_graph(g6, 10.0) == eval_graph(g1, 10.0)

    # Normalization

    degopt8 = normalize!(deepcopy(degopt1))
    (g8, _) = graph_degopt(degopt8)
    (g1, _) = graph_degopt(degopt1)
    @test eval_graph(g8, pi) ≈ eval_graph(g1, pi)

    # get_degopt_coeff
    (HA, HB, y) = get_degopt_coeffs(g8)
    (g9, _) = graph_degopt(Degopt(HA, HB, y))
    @test eval_graph(g8, pi) == eval_graph(g9, pi)
end
