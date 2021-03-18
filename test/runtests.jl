using GraphMatFun
using Test

@testset "GraphMatFun.jl" begin
    include("denman_beavers.jl");
    include("ps_poly.jl");
    include("monomial_poly.jl");
    include("merge_graphs.jl");
    include("graph_ops.jl");
    include("code_gen_test.jl");
    include("bbc.jl");
    include("bbc_sid_exp.jl");
    include("compress.jl");
end
