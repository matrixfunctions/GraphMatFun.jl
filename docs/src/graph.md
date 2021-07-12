
# Graphs and their manipulation


## The computational graph -- Compgraph

The computational graph is represented by a struct `Compgraph{T}`, where `T` is the type of the coefficients.
Each node in the graph is represented by a `Symbol`, and these node-IDs need to be unique for every node.

```@docs
Compgraph
```


## Manipulating graph topology
In general a graph is composed of linear combinations ``C = \alpha A + \beta B`` ([`add_lincomb!`](@ref)), multiplication ``C = AB`` ([`add_mult!`](@ref)), and left inverse ``C = A^{-1}B``  ([`add_ldiv!`](@ref)).
```@docs
add_lincomb!
add_mult!
add_ldiv!
add_sum!
```

The graph can have multiple outputs and the behavior controlled.
```@docs
add_output!
del_output!
clear_outputs!
```

It is also possible to merge two graphs, and to rename and delete nodes.
```@docs
merge_graphs
rename_node!
del_node!
```


## Evaluating the graph
```@docs
eval_graph
get_topo_order
```


## Nodes and coefficients
```@docs
get_sorted_keys
get_all_cref
get_coeffs
set_coeffs!
```


## Errors and error bounds
```@docs
compute_bwd_theta
compute_fwd_theta
compute_bnd_rel_bwd_err
eval_runerr
```


## Reading and writing graphs to file
```@docs
export_compgraph
import_compgraph
```


## Helpers and convenience functions
```@docs
get_polynomial
get_polynomial_coefficients
eltype
big
complex
```
