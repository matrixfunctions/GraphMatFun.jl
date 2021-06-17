using LinearAlgebra
@testset "Linear fit" begin

function graph_test_graph()
    (graph,cref) = graph_monomial(zeros(2))
    add_ldiv!(graph, :invA, :A, :I)
    add_lincomb!(graph, :Gout, 1.0, :P2, 0.0, :invA)
    graph.outputs[1] = :Gout
    push!(cref,(:Gout,2))
    return (graph,cref)
end

    # Polynomial
    coeff = [3; -1; 2.0; 0.1; 1.2]
    p = z -> coeff[1]*one(z) + coeff[2]*z + coeff[3]*z^2 + coeff[4]*z^3 + coeff[5]*z^4
    discr = collect(1:5.0)

    (graph,cref) = graph_monomial(zeros(length(coeff)))
    opt_linear_fit!(graph, p, discr, cref,
                    errtype=:relerr, linlsqr=:svd, droptol=1e-15)
    @test all(abs.(eval_graph(graph,discr) .- p.(discr)) .< 1e-14* abs.(p.(discr)))
    A = [3 4 ; 5 6.6]; PA = p(A)
    @test norm(eval_graph(graph,A) - PA) < 1e-14 * norm(PA)


    (graph,cref) = graph_monomial(zeros(length(coeff)))
    opt_linear_fit!(graph, p, complex(discr), cref,
                    errtype=:relerr, linlsqr=:real_svd, droptol=1e-15)
    @test all(abs.(eval_graph(graph,discr) .- p.(discr)) .< 1e-14* abs.(p.(discr)))
    A = [3 4 ; 5 6.6]; PA = p(A)
    @test norm(eval_graph(graph,A) - PA) < 1e-14 * norm(PA)
    @test isequal(eltype(graph),real(eltype(graph)))
    

    # Other functions
    coeff = [3; -1; 2.0]
    p = z -> coeff[1]*one(z) + coeff[2]*z + coeff[3]*z^(-1)
    discr = [1; -1.1; 1.4]

    (graph,cref) = graph_test_graph()
    opt_linear_fit!(graph, p, discr, cref)
    @test all(abs.(eval_graph(graph,discr) .- p.(discr)) .< 1e-15)


    (graph,cref) = graph_test_graph()
    discr = complex.(discr)
    opt_linear_fit!(graph, p, discr, cref, linlsqr=:real_backslash)
    @test all(abs.(eval_graph(graph,discr) .- p.(discr)) .< 1e-15)
    @test isequal(eltype(graph),real(eltype(graph)))


    # Overdetermined
    (graph,cref) = graph_test_graph()
    discr = collect(1:10.0)
    opt_linear_fit!(graph, p, discr, cref,
                    errtype=:relerr, linlsqr=:nrmeq)
    PA = p(A)
    @test norm(eval_graph(graph,A) - PA) < 1e-13 * norm(PA)


    (graph,cref) = graph_test_graph()
    coeff = complex.(coeff)
    p = z -> coeff[1]*one(z) + coeff[2]*z + coeff[3]*z^(-1)
    opt_linear_fit!(graph, p, complex.(discr), cref,
                    errtype=:relerr, linlsqr=:real_nrmeq)
    @test norm(eval_graph(graph,A) - PA) < 1e-13 * norm(PA)
    @test isequal(eltype(graph),real(eltype(graph)))

end
