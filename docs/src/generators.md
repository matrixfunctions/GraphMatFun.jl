
# Built-in graphs

## Polynomails
Standard evaluation shcemes such as monomial evaluation ([`graph_monomial`](@ref)), Horner evaluation ([`graph_horner`](@ref)) and Paterson--Stockmeyer ([`graph_ps`](@ref)) are available. There are also the degree-optimal polynomials ([`Degopt`](@ref)) described in the next section.
```@docs
graph_monomial
graph_monomial_degopt
graph_horner
graph_horner_degopt
graph_ps
graph_ps_degopt
graph_sastre_poly
```


## Degree-optimal polynomails
The degree-optimal polynomial is a multiplication-economic scheme for evaluating polynomials.
It has the possibility to reach the highest attainable degree possible for a fixed number of multiplications, i.e., degree equal to ``2^m`` with ``m`` multiplication. However, the set of degree-optimal polynomials does not span the whole set of polynomials of degree less than or equal to ``2^m``. Provided [optimization techniques](optim.md) are therefore useful to achieve good approximations.
```@docs
Degopt
graph_degopt
grow!
scale!
square!
normalize!
get_degopt_coeffs
get_degopt_crefs
get_topo_order_degopt
graph_sastre_yks_degopt
```

## Rational functions
```@docs
graph_rational
```

## Matrix exponential
One of the most important matrix functions is the matrix exponential.
The package contains graph-representations for several of the state-of-the-art evaluation schemes.
```@docs
graph_exp_native_jl
graph_exp_native_jl_degopt
graph_sastre_exp
graph_sid_exp
graph_bbc_exp
graph_bbcs_cheb_exp
```

## Other matrix functions and iterations
```@docs
graph_denman_beavers
graph_newton_schulz
graph_newton_schulz_degopt
```