module GraphMatFun

include("compgraph_struct.jl")
include("compgraph_struct_IO.jl")
include("eval.jl")


# Generators
include("generators/denman_beavers_gen.jl");
include("generators/ps_poly.jl");


# Code generations

include("code_gen.jl")

end
