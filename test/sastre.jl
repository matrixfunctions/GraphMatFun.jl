
@testset "Sastre" begin



    for (k,Ak) = [(3,1.1), (4,2.9), (8,17.3)] # Where do I find these limits?
        A = randn(100,100)/20 * Ak
        (graph,cref) = gen_sastre_basic_exp(k)
        @test eval_graph(graph,A) ≈ exp(A)
        @test sum(values(graph.operations) .== :mult) == k
    end
end
