using LinearAlgebra
@testset "compress" begin

    A=[3 4.0; 0.1 -0.3];


    ## Test dangling node removal
    graph=Compgraph();
    add_lincomb!(graph,:T1,1.0,:I,3.0,:A);
    add_mult!(graph,:T2,:T1,:T1);

    # This node is dangling
    add_mult!(graph,:T3,:T2,:A);

    add_output!(graph,:T2);

    graph1=graph;
    graph2=deepcopy(graph);

    compress_graph_dangling!(graph2);
    @test eval_graph(graph1,A) == eval_graph(graph2,A)



    ## Test zero coeff removal (and dangling)
    graph=Compgraph();
    add_mult!(graph,:A2,:A,:A);
    add_mult!(graph,:A4,:A2,:A2);
    add_mult!(graph,:A8,:A4,:A4);
    add_sum!(graph,:T,[3;4;0.1;0;0],[:I;:A;:A2;:A4;:A8]);
    add_output!(graph,:T);

    graph1=graph;
    graph2=deepcopy(graph);
    compress_graph!(graph2);
    @test eval_graph(graph1,A) == eval_graph(graph2,A)

    # Check that we actually removed something
    @test size(get_all_cref(graph2),1) < size(get_all_cref(graph1),1)



    ## Test trivial node detection and removal
    graph=Compgraph();
    add_mult!(graph,:AI,:A,:I);
    add_output!(graph,:AI)
    @test has_trivial_nodes(graph) == true

    graph=Compgraph();
    add_mult!(graph,:IA,:I,:A);
    add_output!(graph,:IA)
    @test has_trivial_nodes(graph) == true

    graph=Compgraph();
    add_ldiv!(graph,:IinvA,:I,:A);
    add_output!(graph,:IinvA)
    @test has_trivial_nodes(graph) == true

    graph=Compgraph();
    add_mult!(graph,:I2,:I,:I); # This node is trivial
    add_mult!(graph,:AI,:A,:I2); # This node is trivial
    add_mult!(graph,:IA,:I2,:AI); # This node is trivial
    add_ldiv!(graph,:P0,:I,:IA); # This node is trivial
    add_output!(graph,:P0);

    graph1=graph;
    graph2=deepcopy(graph);
    @test has_trivial_nodes(graph1) == true
    compress_graph_trivial!(graph2);
    @test eval_graph(graph1,A) == eval_graph(graph2,A)

    # Check that all trivial nodes were removed
    @test isempty(graph2.operations)
    @test isempty(graph2.parents)
    @test isempty(graph2.coeffs)
    @test graph2.outputs == [:A]


end
