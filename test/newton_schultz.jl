using LinearAlgebra
@testset "Newton–Schultz iteration" begin

    A = [2.7 4 ; 5.3 6.6]/7

    (graph,cref) = gen_newton_schulz(4)
    Vk = A
    for i = 1:4
        Vkk = Vk*(2*I - A*Vk)
        Vk = Vkk
    end
    @test eval_graph(graph,A) ≈ Vk

    (graph,cref) = gen_newton_schulz(15)
    @test eval_graph(graph,A) ≈ inv(A)

end
