using GraphMatFun
using Test

@testset "GraphMatFun.jl" begin
    include("denman_beavers.jl");
    include("ps_poly.jl");
    include("graph_ops.jl");
    include("code_gen_test.jl");
    include("bbc.jl");
    include("compress.jl");
end
