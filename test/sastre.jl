
@testset "Sastre" begin



    for (k,Ak) = [(3,1.1), (4,2.9), (6,0.001)] # Where do I find these limits?
        display(k)
        display(Ak)
        A = randn(100,100)/20 * Ak
        (graph,cref) = gen_sastre_basic_exp(k)
        @test eval_graph(graph,A) ≈ exp(A)
        @test sum(values(graph.operations) .== :mult) == k
    end
end
