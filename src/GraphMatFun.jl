module GraphMatFun

include("compgraph_struct.jl")
include("compgraph_struct_IO.jl")
include("eval.jl")
include("compress_graph.jl")
include("merge.jl")
include("degopt.jl")

# Optimization routines
include("optimization/opt_common.jl")
include("optimization/gauss_newton.jl")
include("optimization/linear_fit.jl")

# Generators
include("generators/denman_beavers_gen.jl")
include("generators/ps_poly.jl")
include("generators/degopt_poly.jl")
include("generators/bbc_exp.jl")
include("generators/bbcs_cheb_exp.jl")
include("generators/sid.jl")
include("generators/monomial_poly.jl")
include("generators/rational.jl")
include("generators/horner.jl")
include("generators/newton_schulz.jl")
include("generators/exp.jl")
include("generators/sastre.jl")

# Error bounds
include("error_bounds.jl")

# Code generations
include("code_gen/gen_code.jl")

# Visualization
include("visualization.jl");
end
