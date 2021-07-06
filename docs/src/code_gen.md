
# Code generation and efficiency
While the graph can be evaluated with [`eval_graph`](@ref), the ideal way to utilize a good graph is
to produce efficient code. Currently code generation to [Julia](https://julialang.org/), [Matlab](http://matlab.com), and [C](https://www.iso.org/standard/74528.html) is supported.

## Improving graph topology
The [`Compgraph`](@ref) is a general framework, and sometimes the topology can be improved.
This is especially important before code generation.
```@docs
extract_sums
has_trivial_nodes
has_identity_lincomb
compress_graph! 
compress_graph_output_cleaning!
compress_graph_dangling!
compress_graph_redundant!
compress_graph_trivial!
compress_graph_zero_coeff!
compress_graph_passthrough!
```

## Code generation
```@docs
gen_code
LangJulia
LangMatlab
LangC_MKL
LangC_OpenBLAS
```

