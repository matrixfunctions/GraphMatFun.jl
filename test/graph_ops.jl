using LinearAlgebra, Polynomials
@testset "graphops" begin
    graph = Compgraph()

    # Graph evaluation.

    # Instructions for  (I + 2*(I + A/5))^{-1} (I + 1/2 *(I + A/5)^2)
    add_lincomb!(graph, :C, 1.0, :I, 1 / 5, :A)
    add_mult!(graph, :Csqr, :C, :C)
    add_lincomb!(graph, :N, 1.0, :I, 1 / 2, :Csqr)
    add_lincomb!(graph, :D, [1.0, 2, -2], [:I, :C, :Csqr])
    add_ldiv!(graph, :out0, :D, :N)
    add_lincomb!(graph, :out, [2.0], [:out0]); # Single term
    add_output!(graph, :out)

    A = 0.3
    E = eval_graph(graph, A)
    Z = 2 * ((I + 2 * (I + A / 5)- 2 * (I + A / 5)^2) \ (I + 1 / 2 * (I + A / 5)^2))
    @test Z ≈ E

    A = [0.5; 0.4]
    E = eval_graph(graph, A)
    Z = zero(A)
    for i = 1:2
        a = A[i]
        Z[i] = 2 * ((I + 2 * (I + a / 5)- 2 * (I + a / 5)^2) \ (I + 1 / 2 * (I + a / 5)^2))
    end
    @test Z ≈ E

    A = [0.4 0.2; 1.1 0.3]
    E = eval_graph(graph, A)
    Z = 2 * ((I + 2 * (I + A / 5)- 2 * (I + A / 5)^2) \ (I + 1 / 2 * (I + A / 5)^2))
    @test Z ≈ E

    # Node removal.
    del_node!(graph, :D) #Remove rational part
    del_node!(graph, :out0)
    del_output!(graph, :out)
    @test isempty(graph.outputs)
    del_node!(graph, :out)
    add_output!(graph, :N)
    rename_node!(graph, :N, :N)
    rename_node!(graph, :N, :N_new)
    A = Polynomial("x")
    E = eval_graph(graph, A)
    Z = (1 + 1 / 2 * (1 + A / 5)^2)
    @test Z ≈ E
    del_output!(graph, :N_new)
    @test isempty(graph.outputs)

    # Valid names and manipulation of the output list
    names =
        [:Acoeff_type, :coeff1A, :Bcoeff2, :My_output, :output1_My, :output1MY2]
    for name in names
        add_lincomb!(graph, name, 1.0, :I, 1 / 5, :A)
        add_mult!(graph, name, :I, :A)
        add_ldiv!(graph, name, :D, :N)
        add_output!(graph, name)
    end

    # Clearance of output nodes.
    clear_outputs!(graph)
    @test isempty(graph.outputs)

    # Node renaming
    graph = Compgraph()
    add_lincomb!(graph, :AI, 2.0, :A, 2.0, :I)
    add_mult!(graph, :Pout, :AI, :A)
    add_output!(graph, :Pout)
    try
        rename_node!(graph, :Z, :Z_new)
    catch e
        @test e isa Exception
        @test sprint(showerror, e) == "Node Z not present in the graph."
    end

    graph1 = deepcopy(graph)
    rename_node!(graph1, :AI, :B)
    @test eval_graph(graph1, A) == eval_graph(graph, A)

    graph1 = deepcopy(graph)
    rename_node!(graph1, :Pout, :B)
    @test eval_graph(graph1, A) == eval_graph(graph, A)

    graph1 = deepcopy(graph)
    rename_node!(graph1, :A, :A2tmp)
    add_mult!(graph1, :A2tmp, :A, :A)
    A = [0.4 0.2; 1.1 0.3]
    @test eval_graph(graph1, A) == eval_graph(graph, A^2)

    # Input kwarg to eval_graph.
    graph1 = deepcopy(graph)
    new_input = :B
    rename_node!(graph, :A, new_input)
    @test eval_graph(graph1, A) == eval_graph(graph, A, input = new_input)

    # Removal of output node
    del_output!(graph, :Pout)
    @test isempty(graph.outputs)

    # Setting and getting coefficients
    graph = Compgraph()
    add_lincomb!(graph, :B, [1.0], [:A])
    add_output!(graph, :B)
    set_coeffs!(graph, 2.0, get_all_cref(graph)[1])
    @test graph.coeffs[:B][1] == 2.0
    try
        set_coeffs!(graph, [1.0, 2.0])
    catch e
        @test e isa Exception
        @test sprint(showerror, e) == "Vector of coefficients and defined " *
            "set of coefficients do not have the same length."
    end
    @test get_coeffs(graph, get_all_cref(graph)[1]) == [2.0]

end
