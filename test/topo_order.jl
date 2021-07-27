using LinearAlgebra
@testset "topo order" begin
    g = Compgraph()
    add_lincomb!(g, :X1, 0.2, :I, 0.2, :A)
    add_lincomb!(g, :X2, 0.2, :I, 0.22, :A)
    add_lincomb!(g, :X3, 0.3, :I, 0.1, :A)

    add_mult!(g, :X4, :I, :X3)

    add_mult!(g, :Z1, :X1, :X2)
    add_mult!(g, :Z2, :X2, :X3)
    add_mult!(g, :Z3, :Z1, :Z2)
    add_mult!(g, :Z4, :Z3, :Z3)

    helper = Dict()
    helper[:X1] = -10000 + 1
    helper[:X2] = -10000 + 2
    helper[:X3] = -10000 + 3
    helper[:X4] = -3

    add_output!(g, :Z4)
    order = get_topo_order(g, priohelp = helper, will_not_dealloc = [])[1]
    # Check that 1) helper prioritizes x1,x2,x3 first
    # 2) that x4 comes next since it "frees" :I (which is not in will_not_dealloc)
    @test order[1:4] == [:X1; :X2; :X3; :X4]
    # 3) that z1 comes next since it "frees" if will_not_dealloc=:I
    #    (and before z2 since it has better priority)
    helper[:Z1] = -2
    helper[:Z2] = -1
    order = get_topo_order(g, priohelp = helper, will_not_dealloc = [:I])[1]
    @test order[1:4] == [:X1; :X2; :X3; :Z1]
end
