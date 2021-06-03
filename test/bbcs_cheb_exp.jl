using LinearAlgebra
@testset "BBCS Chebyshev exponential" begin

    for (k,θ) = enumerate([1.38e-5, 2.92e-3, 0.1295, 0.636])#, 2.212, 2*2,212])
        Aθ = θ*I*1im + [0 1e-8; -1e-8 0]
        (graph,cref) = gen_bbcs_cheb_exp(k)
        tv = range(-θ,θ,length=41)
        @test all(abs.(eval_graph(graph,big.(tv))-exp.(-1im*big.(tv))) .<eps()/2 )
        @test eval_graph(graph,Aθ) ≈ exp(-1im*Aθ)
        @test sum(values(graph.operations) .== :mult) == k
    end

end
