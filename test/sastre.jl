
@testset "Sastre" begin

    # Test approximation of exponential and number of multiplications
    for (k,Ak) = [(3,1.1), (4,2.9), (8,17.3)] # Where do I find these limits?
        A = randn(100,100)/20 * Ak
        (graph,cref) = gen_sastre_basic_exp(k)
        @test eval_graph(graph,A) ≈ exp(A)
        @test sum(values(graph.operations) .== :mult) == k
    end


    # Test the degree-9 solver gen_sastre_basic(b) by comparison to tabulated values
    (graph1,cref1) = gen_sastre_basic_exp(3)
    (graph2,cref2) =  gen_sastre_basic(1 ./factorial.(0:8))
    @test sum(values(graph2.operations) .== :mult) == 3
    @test all(cref1 .== cref2)
    for c = cref1
        @test graph1.coeffs[c[1]][c[2]] ≈ graph2.coeffs[c[1]][c[2]]
    end
end
