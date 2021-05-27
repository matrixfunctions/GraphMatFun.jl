using LinearAlgebra, Polynomials

@testset "graphops" begin
    graph = Compgraph();

    # instructions for  (I + 2*(I + A/5))^{-1} (I + 1/2 *(I + A/5)^2)
    add_lincomb!(graph,:C, 1.0, :I, 1/5, :A);
    add_mult!(graph,:Csqr, :C, :C);
    add_lincomb!(graph,:N, 1.0, :I, 1/2, :Csqr);
    add_lincomb!(graph,:D, 1.0, :I, 2, :C);
    add_ldiv!(graph,:out,:D,:N);
    add_output!(graph,:out)

    A = 0.3
    E=eval_graph(graph,A);
    Z = (I+2*(I+A/5))\(I + 1/2 *(I + A/5)^2);
    @test Z≈E

    A = [0.5 ; 0.4];
    E=eval_graph(graph,A);
    Z = zero(A);
    for i=1:2
        a = A[i]
        Z[i] = (I+2*(I+a/5))\(I + 1/2 *(I + a/5)^2);
    end
    @test Z≈E

    A = [0.4 0.2; 1.1 0.3]
    E=eval_graph(graph,A);
    Z = (I+2*(I+A/5))\(I + 1/2 *(I + A/5)^2);
    @test Z≈E

    del_node!(graph,:D) #Remove rational part
    del_output!(graph,:out)
    @test isempty(graph.outputs)
    del_node!(graph,:out)
    add_output!(graph,:N)
    rename_node!(graph,:N,:N)
    rename_node!(graph,:N,:N_new)
    A = Polynomial("x")
    E=eval_graph(graph,A);
    Z = (1 + 1/2 *(1 + A/5)^2);
    @test Z≈E
    del_output!(graph,:N_new)
    @test isempty(graph.outputs)



    # Valid names and manipulation of the output list
    names = [:Acoeff_type, :coeff1A, :Bcoeff2,
             :My_output, :output1_My, :output1MY2]
    for name in names
        add_lincomb!(graph,name, 1.0, :I, 1/5, :A);
        add_mult!(graph,name,:I,:A);
        add_ldiv!(graph,name,:D,:N);
        add_output!(graph,name)
    end
    clear_outputs!(graph)
    @test isempty(graph.outputs)



    # Node renaming
    graph=Compgraph()
    add_lincomb!(graph,:AI,2.0,:A,2.0,:I)
    add_mult!(graph,:Pout,:AI,:A)
    add_output!(graph,:Pout)

    graph1=deepcopy(graph)
    rename_node!(graph1,:AI,:B)
    @test eval_graph(graph1,A)==eval_graph(graph,A)

    graph1=deepcopy(graph)
    rename_node!(graph1,:Pout,:B)
    @test eval_graph(graph1,A)==eval_graph(graph,A)

    graph1=deepcopy(graph)
    rename_node!(graph1,:A,:A2tmp)
    add_mult!(graph1,:A2tmp,:A,:A)
    A = [0.4 0.2; 1.1 0.3]
    @test eval_graph(graph1,A)==eval_graph(graph,A^2)

    del_output!(graph,:Pout)
    @test isempty(graph.outputs)

end
