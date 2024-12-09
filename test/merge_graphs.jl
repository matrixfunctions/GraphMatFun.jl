using LinearAlgebra
@testset "merge" begin
    A = [3 4; 5 6.6]

    # Example which constructs a composition
    x = [3.0; 2.2; 0.1]
    # Generate a polynomial
    (graph1, cref) = graph_ps(x)
    g1_output = graph1.outputs[end]

    # Generate graph x^2
    (graph2, cref) = graph_monomial([0.0; 0; 1.0], input = :B)

    graph = merge_graphs(graph1, graph2, prefix2 = "Monomial")

    # Point the output of graph1 to the input of graph2
    add_lincomb!(graph, :MonomialB, 1.0, g1_output, 0, :I)

    # Adjust outputs
    newout = graph.outputs[end]
    empty!(graph.outputs)
    add_output!(graph, newout)

    @test (x[1] * I + x[2] * A + x[3] * A^2)^2 ≈ eval_graph(graph, A)

    A = [2.7 4; 5.3 6.6] / 7
    (graph, cref) = graph_newton_schulz(3)
    (graph_2, cref_2) = graph_newton_schulz_degopt(3)
    rename_node!(graph_2, graph_2.outputs[1], :Myout)
    (graph_3, cref_3) = graph_newton_schulz_degopt(3, input = :G2Myout)
    graph_4 = merge_graphs(
        graph_2,
        graph_3;
        prefix1 = "G2",
        prefix2 = "G3",
        input2 = :G2Myout,
    )
    @test eval_graph(graph_4, A) ≈ eval_graph(graph, eval_graph(graph, A))
end
@testset "split_lincomb" begin

    A=[
        0.196064  -1.18798    0.0        0.0         0.0         0.0       0.0
        0.342678   0.793844   0.684638   0.0         0.0         0.0       0.0
        1.26005    0.819614  -1.21473   -0.16794     0.0         0.0       0.0
        -0.10478    1.50713    0.138891   0.652377   -0.580547    0.0       0.0
        -0.803933  -1.96458    0.679152  -0.0471995  -0.0454885  -0.715895  0.0
        -0.327724  -1.89847   -0.435685   1.44764    -0.928145    3.06844   1.19071];
    B=[1.46284     0.594921   0.0         0.0        0.0        0.0        0.0
       0.71404    -0.492136  -0.481715    0.0        0.0        0.0        0.0
       0.237897   -0.946491   0.0070323   1.87826    0.0        0.0        0.0
       0.0444275  -1.77078   -0.772886    0.89719    1.55139    0.0        0.0
       -0.310057   -1.13217    1.40868     0.864775  -1.14717   -0.662986   0.0
       -0.160165   -0.525442   0.676407   -0.934459   0.259851  -0.608684  -0.381333];
    c=[-0.5210883864869685
       -0.5011935108466905
       0.9989082416364822
       1.1643621428307982
       -1.4603489516355923
       2.405378746857791
       1.7841340060684932
       0.233340850959019];

    degopt=Degopt(A,B,c);
    n=size(A,1);
    (g,cref)=graph_degopt(degopt);
    g_org=deepcopy(g);
    node=cref[50][1];;
    ind2=[3;6];
    (g,crefs,new_crefs)=split_lincomb!(g,node,ind2;
                                       newnode=Symbol("$(node)_new"),
                                       cref_list=[])

    # Check that it is unmodified
    @test  eval_graph(g,0.1) ≈ eval_graph(g_org,0.1)


    # Check that changing a variable has the same effect
    set_coeffs!(g_org,[3.3],[(node,ind2[1])])
    set_coeffs!(g,[3.3],[new_crefs[1]])

    @test eval_graph(g,0.3) ≈ eval_graph(g_org,0.3)

end
