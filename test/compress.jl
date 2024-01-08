using LinearAlgebra
@testset "compress" begin
    A = [3 4.0; 0.1 -0.3]

    ## Test dangling node removal
    graph = Compgraph()
    add_lincomb!(graph, :T1, 1.0, :I, 3.0, :A)
    add_mult!(graph, :T2, :T1, :T1)

    # This node is dangling
    add_mult!(graph, :T3, :T2, :A)

    add_output!(graph, :T2)

    graph1 = graph
    graph2 = deepcopy(graph)

    compress_graph_dangling!(graph2)
    @test eval_graph(graph1, A) == eval_graph(graph2, A)

    ## Test zero coeff removal (and dangling)
    graph = Compgraph()
    add_mult!(graph, :A2, :A, :A)
    add_mult!(graph, :A4, :A2, :A2)
    add_mult!(graph, :A8, :A4, :A4)
    add_lincomb!(graph, :T, [3; 4; 0.1; 0; 0], [:I; :A; :A2; :A4; :A8])
    add_output!(graph, :T)

    graph1 = graph
    graph2 = deepcopy(graph)
    compress_graph!(graph2)
    @test eval_graph(graph1, A) == eval_graph(graph2, A)

    # Check that we actually removed something
    @test size(get_all_cref(graph2), 1) < size(get_all_cref(graph1), 1)

    ## Test function to detect if graph has linear combination of identities
    graph = Compgraph()
    add_lincomb!(graph, :P1, 1.0, :A, 1.0, :I)
    add_lincomb!(graph, :P2, 1.0, :I, 1.0, :I)
    add_mult!(graph, :P0, :P1, :P2)
    @test has_identity_lincomb(graph) == true

    ## Test trivial node detection and removal
    graph = Compgraph()
    add_mult!(graph, :AI, :A, :I)
    add_output!(graph, :AI)
    @test has_trivial_nodes(graph) == true

    graph = Compgraph()
    add_mult!(graph, :IA, :I, :A)
    add_output!(graph, :IA)
    @test has_trivial_nodes(graph) == true

    graph = Compgraph()
    add_ldiv!(graph, :IinvA, :I, :A)
    add_output!(graph, :IinvA)
    @test has_trivial_nodes(graph) == true

    graph = Compgraph()
    add_mult!(graph, :I2, :I, :I) # This node is trivial
    add_mult!(graph, :AI, :A, :I2) # This node is trivial
    add_mult!(graph, :IA, :I2, :AI) # This node is trivial
    add_ldiv!(graph, :P0, :I, :IA) # This node is trivial
    add_output!(graph, :P0)
    graph1 = graph
    graph2 = deepcopy(graph)
    @test has_trivial_nodes(graph1) == true
    cref = get_all_cref(graph) # Empty as graph has no linear combination nodes.
    compress_graph_trivial!(graph2, cref)
    @test isempty(cref)
    @test eval_graph(graph1, A) == eval_graph(graph2, A)
    @test isempty(graph2.operations)
    @test isempty(graph2.parents)
    @test isempty(graph2.coeffs)
    @test graph2.outputs == [:A]

    graph = Compgraph()
    add_lincomb!(graph, :P1, [1.0,2.0,0], [:A, :I,:A])
    add_mult!(graph, :P2a, :A, :I)
    add_mult!(graph, :P2b, :I, :A)
    add_mult!(graph, :P2, :P2a, :P2b)
    add_mult!(graph, :P0, :P1, :P2)
    #add_lincomb!(graph, :Px,0,:I,1.0,:P0)
    add_output!(graph, :P0)

    graph1 = graph
    graph2 = deepcopy(graph)
    @test has_trivial_nodes(graph1) == true
    cref = get_all_cref(graph) # Empty as graph has no linear combination nodes.
    compress_graph_trivial!(graph2, cref)
    @test graph2.outputs == [:P0]

    ## Test function to remove redundant nodes.
    graph = Compgraph()
    add_mult!(graph,:A2,:A,:A);
    add_lincomb!(graph, :AI1, [1.0,2.0,3.0], [:A, :I, :A2])
    add_lincomb!(graph, :AI2, [1.0,2.0,3.0], [:A, :I, :A2])
    add_lincomb!(graph, :IA1, 2.0, :I, 1.0, :A)
    add_mult!(graph, :P1, :AI1, :AI2)
    add_mult!(graph, :P2, :AI2, :AI1)
    add_ldiv!(graph, :D1, :AI1, :P2)
    add_ldiv!(graph, :D2, :AI1, :P2)
    add_mult!(graph, :out, :D1, :D2)
    add_output!(graph, :out)

    cref = get_all_cref(graph)

    # Reduce redundant multiplications and divisions
    graph1 = deepcopy(graph)
    cref1 = get_all_cref(graph1)
    compress_graph_redundant!(graph1, cref, compress_lincomb = false)
    @test length(get_sorted_keys(graph1)) == 7
    @test cref == cref1
    @test eval_graph(graph1, A) == eval_graph(graph, A)

    # Reduce redundant multiplications and lincombs
    graph1 = deepcopy(graph)
    cref1 = get_all_cref(graph1)
    compress_graph_redundant!(graph1, cref1, compress_lincomb = true)
    @test length(get_sorted_keys(graph1)) == 6
    @test Set(get_all_cref(graph1))==Set(cref1)
    @test eval_graph(graph1, A) == eval_graph(graph, A)



    # Zero coeff removal and passthrough
    graph = Compgraph()
    add_mult!(graph,:A2,:A,:A);
    add_lincomb!(graph, :Ax, [0,0.0,1.0], [:A, :I, :A2])
    add_ldiv!(graph, :out,:Ax,:A);
    add_output!(graph,:out);

    cref=get_all_cref(graph);
    graph1=deepcopy(graph);
    compress_graph_zero_coeff!(graph1,cref)
    compress_graph_passthrough!(graph1,cref)
    @test eval_graph(graph,A) == eval_graph(graph1,A)

end
