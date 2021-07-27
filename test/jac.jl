@testset "Jacobian" begin
    graph = Compgraph()
    c = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    add_lincomb!(graph, :B1, c[1], :I, c[2], :A)
    add_lincomb!(graph, :B2, c[3], :I, c[4], :A)
    add_mult!(graph, :N, :B1, :B2)
    add_lincomb!(graph, :D, c[5], :I, c[6], :A)
    add_ldiv!(graph, :R, :D, :N)
    add_output!(graph, :R)

    x = randn(5)
    @test eval_graph(graph, x) ≈
          (c[1] .+ c[2] * x) .* (c[3] .+ c[4] * x) ./ (c[5] .+ c[6] * x)

    dg_dc3 = eval_jac(graph, x, (:B2, 1))
    @test dg_dc3 ≈ (c[1] .+ c[2] * x) ./ (c[5] .+ c[6] * x)

    dg_dc6 = eval_jac(graph, x, (:D, 2))
    @test dg_dc6 ≈
          -x .* (c[1] .+ c[2] * x) .* (c[3] .+ c[4] * x) ./
          (c[5] .+ c[6] * x) .^ 2

    dg_dc2c5 = eval_jac(graph, x, [(:B1, 2), (:D, 1)])
    @test dg_dc2c5 ≈ hcat(
        x .* (c[3] .+ c[4] * x) ./ (c[5] .+ c[6] * x),
        -(c[1] .+ c[2] * x) .* (c[3] .+ c[4] * x) ./ (c[5] .+ c[6] * x) .^ 2,
    )
end
