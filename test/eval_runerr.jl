using LinearAlgebra
@testset "error analyis" begin
    # Eval runerr
    # Same as Example 17 in manuscript
    graph = Compgraph(Float64);
    add_lincomb!(graph, :y, 1, :I, 1, :x)
    add_lincomb!(graph, :z, 1, :I, 1/2, :x)
    add_mult!(graph, :y2, :y, :y)
    add_mult!(graph, :z2, :z, :z)
    add_lincomb!(graph, :out, 1, :y2, -1, :z2)
    add_output!(graph, :out)
    x = 2. .^ -27
    running_error_bound = eval_runerr(graph, x, input = :x)
    @test running_error_bound < 3e-7;

    x = [1., 2.] .^ -27
    running_error_bound = eval_runerr(graph, x, input = :x)
    @test all(running_error_bound .< [2e-7, 3e-7])

    # compute theta
    (graph,_)=graph_ps([1.0,1.0, 1 / 2, 1 /6])
    q=compute_bnd_rel_bwd_err(:exp,graph,numterms=10,numdigits=30)
    @test q(0.1) < 1e-04


end
