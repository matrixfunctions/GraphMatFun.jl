using GraphMatFun
using Test

@testset "GraphMatFun.jl" begin
    include("bbc_sid_exp.jl");
    include("bbc.jl");
    include("bbc_cheb_exp.jl");
    include("code_gen_test.jl");
    include("compress.jl");
    include("degopt.jl");
    include("denman_beavers.jl");
    include("exp.jl");
    include("gauss_newton.jl");
    include("graph_ops.jl");
    include("import_export.jl");
    include("linear_fit.jl");
    include("merge_graphs.jl");
    include("newton_schulz.jl");
    include("opt_common.jl");
    include("polynomials_degopt.jl");
    include("polynomials.jl");
    include("rational.jl");
    include("topo_order.jl");
end
