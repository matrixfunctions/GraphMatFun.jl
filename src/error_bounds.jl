using Polynomials, Roots
export get_polynomial
export get_polynomial_coefficients
export compute_fwd_theta
export compute_bnd_rel_bwd_err
export compute_bwd_theta

"""
    p = get_polynomial(graph::Compgraph)

Return the polynomial underlying the computational `graph`. The `graph` is
assumed to involve only linear combinations and multiplications.

`p` is polynomial of type `Polynomial` from the package
[`Polynomials.jl`](https://juliamath.github.io/Polynomials.jl/stable/).

See also [`get_polynomial_coefficients`](@ref).
"""
function get_polynomial(graph::Compgraph)
    graph = Compgraph(Any, graph) # Allow for symbolic values
    p = eval_graph(graph, Polynomial("x"))
    return p
end

"""
    c = get_polynomial_coefficients(graph::Compgraph)

Return the coefficients of the polynomial underlying the computational `graph`.
The `graph` is assumed to involve only linear combinations and multiplications.

The coefficients `c` is a vector containing the coefficients as expressed in
the monomial basis and sorted from the leading coefficient to the constant term.

See also [`get_polynomial`](@ref).
"""
function get_polynomial_coefficients(graph::Compgraph)
    return get_polynomial(graph).coeffs
end

"""
    (e_fwd,theta)=compute_fwd_theta(
        graph::Compgraph{T},
        f;
        coefftype = T,
        tolerance = eps(coefftype) / 2,
        theta_init = big"0.1",
    )

Return relative forward error and corresponding theta for approximation to `f`.

The first output ``e_{fwd}`` is the function that computes the relative forward
error``e_{fwd}(z)=|f(z) - p(z)| / |z|``, where ``p`` is the polynomial
underlying the computational graph `graph`. This is meaningful only if ``p``
approximates ``f`` in some sense.

The second output ``θ`` is the largest positive real number such that
``e_{fwd}(θ) ≤ tolerance``, where ``tolerance`` is by default the unit roundoff
of the data type `coefftype`, which in turn defaults to the type of the
coefficients of `graph`.

The value of ``θ`` is approximated by using the built-in function `fzero` with
starting value set to `theta_init`. By default, this root-finding procedure
uses high-precision arithmetic.
"""
function compute_fwd_theta(
    graph::Compgraph{T},
    f;
    coefftype = T,
    tolerance = eps(coefftype) / 2,
    theta_init = big"0.1",
) where {T}
    # Obtain coefficients of polynomial.
    coeff = get_polynomial_coefficients(graph)

    # Convert coefficients to required type.
    coeff = convert.(coefftype, coeff)
    if isreal(coeff)
        coeff = real(coeff)
    end

    # Find point where the relative forward error equals tolerance.
    p = Polynomial(abs.(coeff), :x)
    e_fwd(z) = abs.(f.(z) - p.(z)) ./ abs.(z)
    g(z) = e_fwd(z) - tolerance
    theta_fwd = fzero(g, theta_init, verbose = true)
    return e_fwd, convert(coefftype, theta_fwd)
end

"""
    e_bwd=compute_bnd_rel_bwd_err(
        f,
        graph::Compgraph{T};
        coefftype = T,
        numterms = 100,
        numdigits = 100,
    )

Compute a bound on the relative backward error of the function `f`.

The returned function ``e_{bwd}`` is a bound on the relative backward error of
the polynomial underlying the computational graph `graph`, seen as an
approximant to `f`. For the underlying polynomial ``p``, the
function returns a bound on `|δ|` where `δ` is such that ``f(z+δ) = p(z)``.

The first argument `f` is a symbol. Currently supported values include:

  - `:exp` The bound is computed by means of the identity δ = log((exp(-z)
    p(z)-1)+1). The code computes the series expansion of the right-hand side of the
    equation, truncates it to the first `numterms` coefficients, bounds each
    coefficient of the ensuing polynomial with its absolute value. The computation
    uses `numdigits` digits of precision.
"""
function compute_bnd_rel_bwd_err(
    f,
    graph::Compgraph{T};
    coefftype = T,
    numterms = 100,
    numdigits = 100,
) where {T}
    if f != :exp
        error("Currently :exp is the only supported function.")
    end
    # Set precision to numdigits decimal digits.
    setprecision(Integer(ceil(log2(10.0^numdigits))))

    # Obtain coefficients of polynomial p.
    coeff = get_polynomial_coefficients(graph)

    # Convert coefficients of p to required type.
    coeff = convert.(coefftype, coeff)
    if isreal(coeff)
        coeff = real(coeff)
    end

    # Series expansion of polynomial approximant p.
    papproximant = Polynomial(coeff)

    # (Truncated) series expansions of exp(-z) and log(1+z).
    expminusz =
        (-1) .^ (0:1:numterms) ./
        factorial.(collect(big(0.0):big(1.0):numterms))
    pexpminusz = Polynomial(expminusz)
    logzplusone = [0; (-1) .^ (0:1:numterms-1) ./ (1:1:numterms)]
    plogzplusone = Polynomial(logzplusone)

    # From exp(z+δ) ≈ p(z), approximate δ ≈ log((exp(-z) p(z)-1)+1).
    # The variables are as follows:
    #    * pexpzpz: the coefficients of exp(-z) p(z)
    #    * presult: the coefficients of log((exp(-z) p(z)-1)+1)
    pexpzpz = Polynomial((pexpminusz*papproximant).coeffs[1:numterms])
    presult = Polynomial((plogzplusone(pexpzpz - 1)).coeffs[1:numterms])

    # Compute bound on |δ|.
    bnd_bwd_err = Polynomial(abs.(presult.coeffs))
    e_bwd(z) = abs.(bnd_bwd_err(z)) ./ abs.(z)

    return e_bwd
end

"""
    (e_bwd,theta)=compute_bwd_theta(;
        bnd_rel_err = compute_bnd_rel_bwd_err(:exp, graph),
        tolerance = eps(coefftype) / 2,
        theta_init = big"0.2",
        use_log = false,
    )

    (e_bwd,theta)=compute_bwd_theta(
        graph::Compgraph{T};
        coefftype = T,
        numterms = 100,
        numdigits = 100,
        tolerance = eps(coefftype) / 2,
        theta_init = big"0.2",
        use_log = false,
    )

Compute bound on relative backward error with corresponding theta.

The first output ``e_{bwd}`` is a function that returns a bound on the relative
backward error of the polynomial underlying the computational graph `graph`,
seen as an approximant to the exponential. For the underlying polynomial ``p``,
the function returns a bound on `|δ|` where `δ` is such that exp(z+δ) = p(z).
This function is computed using `compute_bnd_rel_bwd_err(:graph)`.

The second output ``θ`` is the largest positive real number such that
``e_{bwd}(θ) ≤ tolerance``, where ``tolerance`` is by default the unit roundoff
of the data type `coefftype`, which in turn defaults to the type of the
coefficients of `graph`.

The value of ``theta`` is estimated by approximately solving the equation
``e_{bwd}(z) = tolerance``, using the built-in function `fzero` with starting
value set to `theta_init`. By default this root-finding procedure uses
high-precision arithmetic.

If the kwarg `use_log` is set to `true`, then the value of ``theta`` is computed
by approximating a solution to the equation ``log(e_{bwd}(z)) = log(tolerance)``
instead.

The alternative form accepts a function that returns a bound on the relative
backward error. By default, the function is constructed with
`compute_bnd_rel_bwd_err(:exp, ...)` and the default values for the kwargs in
the first form.
"""
function compute_bwd_theta(;
    bnd_rel_err = compute_bnd_rel_bwd_err(:exp, graph),
    tolerance = eps(coefftype) / 2,
    theta_init = big"0.2",
    use_log = false,
)
    # Find point where bound on relative backward error equals tolerance.
    e_bwd = bnd_rel_err
    h(z) = use_log ? log(e_bwd(z) / tolerance) : e_bwd(z) - tolerance
    theta_bwd = fzero(h, theta_init)

    return e_bwd, theta_bwd
end

function compute_bwd_theta(
    graph::Compgraph{T};
    coefftype = T,
    numterms = 100,
    numdigits = 100,
    tolerance = eps(coefftype) / 2,
    theta_init = big"0.2",
    use_log = false,
) where {T}
    return compute_bwd_theta(
        bnd_rel_err = compute_bnd_rel_bwd_err(
            :exp,
            graph;
            coefftype = coefftype,
            numterms = numterms,
            numdigits = numdigits,
        ),
        tolerance = tolerance,
        theta_init = theta_init,
        use_log = use_log,
    )
end
