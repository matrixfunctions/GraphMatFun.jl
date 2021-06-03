using LinearAlgebra
@testset "BBCS Chebyshev exponential" begin

    # Basic implementation
    for (k,θ) = enumerate([1.38e-5, 2.92e-3, 0.1295, 0.636, 2.212])
        Aθ = θ*I*1im + [0 1e-8; -1e-8 0]
        (graph,cref) = gen_bbcs_cheb_exp(k,T=Complex{BigFloat})
        tv = range(-θ,θ,length=41)
        @test all(abs.(eval_graph(graph,big.(tv))-exp.(-1im*big.(tv))) .<eps()/2 )
        (graph,cref) = gen_bbcs_cheb_exp(k,T=ComplexF64)
        @test eval_graph(graph,Aθ) ≈ exp(-1im*Aθ)
        @test sum(values(graph.operations) .== :mult) == k
    end

    # Scalign-and-squaring phase
    for (k,θ) = enumerate([2*2.212, 4*2.212, 8*2.212, 16*2.212, 32*2.212])
        kk = k+5
        Aθ = θ*I*1im + [0 1e-8; -1e-8 0]
        (graph,cref) = gen_bbcs_cheb_exp(kk,T=Complex{BigFloat})
        tv = range(-θ,θ,length=41)
        @test all(abs.(eval_graph(graph,big.(tv))-exp.(-1im*big.(tv))) .<eps()*2^(k-1) )
        (graph,cref) = gen_bbcs_cheb_exp(kk,T=ComplexF64)
        @test eval_graph(graph,Aθ) ≈ exp(-1im*Aθ)
        @test sum(values(graph.operations) .== :mult) == kk
    end

end
