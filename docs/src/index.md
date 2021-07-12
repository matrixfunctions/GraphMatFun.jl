```@meta
CurrentModule = GraphMatFun
```

# GraphMatFun.jl

[GraphMatFun.jl](https://github.com/matrixfunctions/GraphMatFun.jl) is a [Julia](https://julialang.org/) package
for working with computational [graphs](https://en.wikipedia.org/wiki/Graph_(abstract_data_type)) representing
[functions of matrices](https://en.wikipedia.org/wiki/Analytic_function_of_a_matrix).


## Installation

The package can be installed with the command
```julia-repl
julia> ]
(v1.7) pkg> add https://github.com/matrixfunctions/GraphMatFun.jl.git
```
When installed, the package can be loaded into the current session by writing
```julia-repl
julia> using GraphMatFun
```


## Building a first graph

Consider the following polynomial,
```math
p(x) = (1+2x+x^2)^2.
```
A graph representing ``p(x)`` can be manually constructed by defining the operations.
```julia-repl
julia> graph = Compgraph(Float64);
julia> add_mult!(graph,:A2,:A,:A);
julia> add_lincomb!(graph,:P1,1.0,:I,2.0,:A);
julia> add_lincomb!(graph,:P2,1.0,:A2,1.0,:P1);
julia> add_mult!(graph,:P,:P2,:P2);
julia> add_output!(graph,:P);
```

The graph can be evaluated in both scalars, vectors, and matrices. For vectors the evaluation is elementwise.
```julia-repl
julia> x = [1, 2, 3.5];
julia> eval_graph(graph, x)
3-element Vector{Float64}:
  16.0
  81.0
 410.0625
julia> (1 .+ 2*x + x.^2).^2
3-element Vector{Float64}:
  16.0
  81.0
 410.0625
julia> A = [1 3; 7 2.0];
julia> eval_graph(graph, A)
2×2 Matrix{Float64}:
 1150.0   825.0
 1925.0  1425.0
julia> (I + 2*A + A^2)^2
2×2 Matrix{Float64}:
 1150.0   825.0
 1925.0  1425.0
```

Since the graph is a polynomial the coefficients can be directly extracted.
```julia-repl
julia> get_polynomial_coefficients(graph)
5-element Vector{Float64}:
 1.0
 4.0
 6.0
 4.0
 1.0
```


## How do I cite it?

We are finalizing a software paper for this work.
If you find this software useful please cite it by using this citation data:
```bibtex
@Misc{,
  author = {E. Jarlebring and M. Fasi and E. Ringh},
  title  = {Computational graphs for matrix functions},
  year   = {2021},
}
```
