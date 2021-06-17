using LinearAlgebra
@testset "Gauss--Newton" begin

    coeff = [3; -1; 2.0; 0.1; 1.2]
    p = z -> coeff[1]*one(z) + coeff[2]*z + coeff[3]*z^2 + coeff[4]*z^3 + coeff[5]*z^4
    discr = collect(5*exp.(1im*range(0,2*pi,length=21))[1:end-1])

    coeff_pert = [3; -1; 2.0; 0.1; 1.2] + 0.9*[-3; 1; -2.0; 0.1; -1.2]
    (graph,cref) = graph_monomial(coeff_pert);
    graph = complex(graph)
    opt_gauss_newton!(graph, p, discr, logger=0,
                      stoptol=1e-14, γ0=0.85 , cref=cref,
                      errtype=:relerr,
                      linlsqr=:svd, droptol=1e-15)

    @test all(abs.(eval_graph(graph,discr) .- p.(discr))
              .< 1e-14* abs.(p.(discr)))
    A = [3 4 ; 5 6.6]; PA = p(A)
    @test norm(eval_graph(graph,A) - PA) < 1e-14 * norm(PA)

    new_input = :B
    (graph,cref) = graph_monomial(coeff_pert, input=new_input);
    graph = complex(graph)
    opt_gauss_newton!(graph, p, discr, logger=0,
                      stoptol=1e-14, input=new_input, γ0=0.85 , cref=cref,
                      errtype=:relerr,
                      linlsqr=:svd, droptol=1e-15)

    @test all(abs.(eval_graph(graph,discr,input=:B) .- p.(discr))
              .< 1e-14* abs.(p.(discr)))
    A = [3 4 ; 5 6.6]; PA = p(A)
    @test norm(eval_graph(graph,A,input=:B) - PA) < 1e-14 * norm(PA)



    coeff_pert = [3; -1; 2.0; 0.1; 1.2] + 0.9*[-3; 1; -2.0; 0.1; -1.2]
    (graph,cref) = graph_monomial(coeff_pert);
    graph = complex(graph)
    opt_gauss_newton!(graph, p, discr,
                      stoptol=5e-13, cref=cref, errtype=:abserr)
    @test all(abs.(eval_graph(graph,discr) .- p.(discr)) .< 5e-13)

    new_input = :Anew
    rename_node!(graph,:A,new_input,cref)
    graph = complex(graph)
    opt_gauss_newton!(graph, p, discr,
                      stoptol=5e-13, input=new_input,
                      cref=cref, errtype=:abserr)
    @test all(abs.(eval_graph(graph,discr,input=new_input) .- p.(discr))
              .< 5e-13)

end
