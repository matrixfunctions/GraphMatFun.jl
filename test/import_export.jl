using LinearAlgebra
@testset "Graph export and import" begin
    fnames = Vector{String}()
    tests = Vector{Any}(undef, 10)
    tests[1:5] = randn(ComplexF64, 5)
    for i = 1:5
        tests[5+i] = randn(ComplexF64, 20, 20)
    end

    den_coeffs = 10 * randn(7)
    num_coeffs = 10 * randn(7)

    (graph, cref) = graph_rational(den_coeffs, num_coeffs, graph_ps)
    add_lincomb!(graph,:Qnew,[1.0, -1.0, 3.3], [graph.outputs[1],:I, :A])
    clear_outputs!(graph);
    add_output!(graph,:Qnew)
    fname = string(tempname(), ".cgr")
    export_compgraph(graph, fname)
    graph2 = import_compgraph(fname)
    for i = 1:10
        @test eval_graph(graph, tests[i]) ≈ eval_graph(graph2, tests[i])
    end
    rm(fname)

    graph = Compgraph(ComplexF64, graph)
    x = get_coeffs(graph)
    x = x + 0.1 * randn(ComplexF64, length(x))
    set_coeffs!(graph, x)
    graph1 = deepcopy(graph)

    add_output!(graph1, cref[1][1])
    add_output!(graph1, cref[6][1])
    add_output!(graph1, cref[12][1])

    fname = string(tempname(), ".cgr")
    order = get_topo_order(graph)[1]
    export_compgraph(
        graph1,
        fname,
        fun = "TEST",
        dom = "[-2,0.1]",
        err = "???",
        genby = basename(@__FILE__),
        order = order,
    )
    graph2 = import_compgraph(fname)
    rm(fname)
    for i = 1:10
        @test eval_graph(graph, tests[i]) != eval_graph(graph1, tests[i])
        @test eval_graph(graph, tests[i]) ≈ eval_graph(graph2, tests[i])
    end
end
