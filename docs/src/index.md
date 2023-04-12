```@meta
CurrentModule = GraphMatFun
```

# GraphMatFun.jl

[GraphMatFun.jl](https://github.com/matrixfunctions/GraphMatFun.jl) is a [Julia](https://julialang.org/) package
for working with computational [graphs](https://en.wikipedia.org/wiki/Graph_(abstract_data_type)) representing
[functions of matrices](https://en.wikipedia.org/wiki/Analytic_function_of_a_matrix).


## Installation

The package is registered with the Julia central registry and can be installed with the command
```julia-repl
julia> ]
(v1.7) pkg> add GraphMatFun
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

If you find this software useful, please cite the open access [article](https://doi.org/10.1145/3568991):
```bibtex
@Article{jfr22,
  author = "Jarlebring, Elias and Fasi, Massimiliano and Ringh, Emil",
  title = "Computational Graphs for Matrix Functions",
  journal = "ACM Trans. Math. Software",
  year = 2023,
  volume = 48,
  number = 4,
  pages = "1-35",
  month = mar,
  doi = "10.1145/3568991"
  }
```
Previous versions of the manuscript are available as [arXiv:2107.12198 [math.NA]](https://arxiv.org/abs/2107.12198) and [MIMS Eprint 2021.12](https://eprints.maths.manchester.ac.uk/2858/).
