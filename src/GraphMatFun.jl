module GraphMatFun

include("compgraph_struct.jl")
include("compgraph_struct_IO.jl")
include("eval.jl")
include("compress_graph.jl")
include("merge.jl")

# Optimization routines
include("optimization/opt_common.jl")
include("optimization/gaussnewton.jl")

# Generators
include("generators/denman_beavers_gen.jl");
include("generators/ps_poly.jl");
include("generators/bbc.jl");
include("generators/bbc_exp.jl");
include("generators/sid.jl");
include("generators/monomial_poly.jl");


# Code generations
include("gen_code.jl")

end
