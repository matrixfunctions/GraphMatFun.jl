using LinearAlgebra
@testset "Matrix exponential" begin


    exp_gens = Dict("Julia native"        => :gen_exp_native_jl,
                    "Julia native degopt" => :gen_exp_native_jl_degopt)#,
                    # "expm 2009"           => :gen_expm2009)

    for key = keys(exp_gens)
        @testset "$key" begin
            for (i,k) = enumerate([0.014, 0.24, 0.94, 2.0, 5.3, 10.7, 21.5, 43.1, 86.0, 172.0, 345.0, 691.0])
                Ak = k*I + [0 1e-8; 0 0]
                (graph,cref) = eval(exp_gens[key])(Ak)
                @test eval_graph(graph,Ak) ≈ exp(Ak)
                if key != "expm 2009"
                    @test sum(values(graph.operations) .== :mult) == i+1
                end

            end
        end
    end
end
