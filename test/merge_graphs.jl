using LinearAlgebra
@testset "merge" begin
    A=[3 4 ; 5 6.6];

    # Example which constructs a composition
    x=[3.0;2.2;0.1];
    # Generate a polynomial
    (graph1,cref)=gen_ps(x);
    g1_output=graph1.outputs[end]

    # Generate graph x^2
    (graph2,cref)=gen_monomial([0.0;0;1.0],input=:B);


    graph=merge_graphs(graph1,graph2,prefix2="Monomial");

    # Point the output of graph1 to the input of graph2
    add_lincomb!(graph,:MonomialB,1.0,g1_output,0,:I);

    # Adjust outputs
    newout=graph.outputs[end];
    empty!(graph.outputs);
    add_output!(graph,newout);

    @test (x[1]*I+x[2]*A+x[3]*A^2)^2≈eval_graph(graph,A)


end
