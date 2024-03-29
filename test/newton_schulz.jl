using LinearAlgebra
@testset "Newton–Schulz iteration" begin
    A = [2.7 4; 5.3 6.6] / 7

    k = 4
    (graph, cref) = graph_newton_schulz(k)
    (graph_2, cref_2) = graph_newton_schulz_degopt(k)
    Vk = A
    for i = 1:k
        Vkk = Vk * (2 * I - A * Vk)
        Vk = Vkk
    end
    @test eval_graph(graph, A) ≈ Vk
    @test eval_graph(graph_2, A) ≈ Vk

    (graph, cref) = graph_newton_schulz(15)
    @test eval_graph(graph, A) ≈ inv(A)

    (graph_2, cref_2) = graph_newton_schulz_degopt(15)
    @test eval_graph(graph, A) ≈ eval_graph(graph_2, A)
end
