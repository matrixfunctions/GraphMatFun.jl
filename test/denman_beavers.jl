@testset "denman_beavers" begin
    A = [0.1 0.3; 0.2 0.7]
    (graph, cref) = graph_denman_beavers(10)
    @test sqrt(A) ≈ eval_graph(graph, A)
end
