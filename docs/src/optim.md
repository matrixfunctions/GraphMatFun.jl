
# Graph optimization
For a fixed topology of the graph, the coefficients can be optimized.
The package provides some routines for doing the optimization.
However, by using the functions [`eval_graph`](@ref) and [`eval_jac`](@ref), external routines based on
function values and first derivatives can in principle be used.

```@docs
opt_gauss_newton!
eval_jac
opt_linear_fit!
solve_linlsqr!
adjust_for_errtype!
```