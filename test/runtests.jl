using GraphMatFun
using Test

@testset "GraphMatFun.jl" begin
    include("bbc.jl");
    include("bbc_sid_exp.jl");
    include("code_gen_test.jl");
    include("compress.jl");
    include("denman_beavers.jl");
    include("gauss_newton.jl");
    include("graph_ops.jl");
    include("linear_fit.jl");
    include("merge_graphs.jl");
    include("newton_schulz.jl");
    include("opt_common");
    include("polynomials.jl");
    include("polynomials_degopt.jl");
end
