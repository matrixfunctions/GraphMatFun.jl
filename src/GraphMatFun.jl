module GraphMatFun

include("compgraph_struct.jl")
include("compgraph_struct_IO.jl")
include("eval.jl")
include("compress_graph.jl")


# Generators
include("generators/denman_beavers_gen.jl");
include("generators/ps_poly.jl");
include("generators/bbc.jl");


# Code generations

include("gen_code.jl")

end
