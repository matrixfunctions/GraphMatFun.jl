
# Graphs and their manipulation
```@docs
Compgraph
```

## Manipulating graph topology
```@docs
add_mult!
add_lincomb!
add_ldiv!
add_sum!
add_output!
del_output!
clear_outputs!
rename_node!
del_node!
merge_graphs
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
compute_fwd_theta
compute_bwd_theta_exponential
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

