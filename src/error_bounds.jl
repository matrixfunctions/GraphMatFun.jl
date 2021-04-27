using Polynomials,Roots
export get_polynomial
export get_polynomial_coefficients
export compute_fwd_theta
export compute_bwd_theta_exponential

"""
    get_polynomial(graph::Compgraph)

Return the polynomial underlying the computational graph.
"""
function get_polynomial(graph::Compgraph)
    graph=Compgraph(Any, graph); # Allow for symbolic values
    p=eval_graph(graph, Polynomial("x"));
    return p
end

"""
    get_polynomial_coefficients(graph::Compgraph)

Return the coefficients of the polynomial underlying the computational graph.
The coefficients are expressed in the monomial basis and sorted from the
leading coefficient to the constant term.
"""
function get_polynomial_coefficients(graph::Compgraph)
    return get_polynomial(graph).coeffs
end

"""
    (e_fwd,theta)=compute_fwd_theta(graph::Compgraph{T},f;
                                    coefftype=T,
                                    tolerance=eps(coefftype)/2
                                    theta_init=big"0.1") where T

Return relative forward error and corresponding theta for approximation to `f`.

The first output ``e_{fwd}`` is the function that computes the relative forward
error``e_{fwd}(z)=|f(z) - p(z)| / |z|``, where ``p`` is the polynomial
underlying the computational graph `graph`. This is meaningful only if ``p`` is
approximates ``f`` in some sense.

The second output ``θ`` is the largest positive real number such that
 ``e_{fwd}(θ) ≤ tolerance``, where ``tolerance`` is by default the unit
 roundoff of the data type `coefftype`, which in turn defaults to the type of
 the coefficients of `graph`.

 The value of ``θ`` is approximated by using the built-in function `fzero` with
 starting value set to `theta_init`. By default, this root-finding procedure
 uses high-precision arithmetic.
"""
function compute_fwd_theta(graph::Compgraph{T},f;
                           coefftype=T,
                           tolerance=eps(coefftype)/2,
                           theta_init=big"0.1") where T
    # Obtain coefficients of polynomial.
    coeff=get_polynomial_coefficients(graph)

    # Convert coefficients to required type.
    coeff=convert.(coefftype,coeff);
    if isreal(coeff)
        coeff=real(coeff)
    end

    # Find point where the relative forward error equals tolerance.
    p=Polynomial(abs.(coeff),:x)
    e_fwd(z)=abs.(f.(z)-p.(z)) ./ abs.(z)
    g(z)=e_fwd(z)-tolerance
    theta_fwd=fzero(g,theta_init,verbose=true)
    return e_fwd,convert(coefftype,theta_fwd)
end

"""
    (e_bwd,theta)=compute_bwd_theta_exponential(graph::Compgraph{T},
                                                nterms=100,
                                                ndigits=100,
                                                coefftype=T,
                                                tolerance=eps(coefftype)/2,
                                                theta_init=big"0.1") where T

Bound on relative backward error of the exponential with corresponding theta.

The first output ``e_{bwd}`` is the function that computes a bound on the
relative backward error of the polynomial underlying the computational graph
`graph`, seen as an approximant to the exponential. In other words, for the
underlying polynomial ``p`` function returns a bound on `|δ|` where `δ` is such
that exp(z+δ) = p(z).

The bound is computed by means of the identity δ = log((exp(-z) p(z)-1)+1). The
code computes the series expansion of the right-hand side of the equation,
truncates it to the first `nterms` coefficients, bounds each coefficient of the
ensuing polynomial with its absolute value. The computation uses `ndigits`
digits of precision.

The second output ``θ`` is the largest positive real number such that
 ``e_{bwd}(θ) ≤ tolerance``, where ``tolerance`` is by default the unit
 roundoff of the data type `coefftype`, which in turn defaults to the type of
 the coefficients of `graph`.

 The value of ``theta`` is approximated by using the built-in function `fzero`
 with starting value set to `theta_init`. By default this root-finding procedure
 uses high-precision arithmetic.
"""
function compute_bwd_theta_exponential(graph::Compgraph{T};
                                       coefftype=T,
                                       nterms=100,
                                       ndigits=100,
                                       tolerance=eps(coefftype)/2,
                                       theta_init=big"0.2") where T
    # Set precision to ndigits decimal digits.
    setprecision(Integer(ceil(log2(10. ^ndigits))));

    # Obtain coefficients of polynomial p.
    coeff=get_polynomial_coefficients(graph)

    # Convert coefficients of p to required type.
    coeff=convert.(coefftype,coeff);
    if isreal(coeff)
        coeff=real(coeff)
    end

    # Series expansion of polynomial approximant p.
    a=Polynomial(coeff)

    # (Truncated) series expansions of exp(-z) and log(1+z).
    expminusz=(-1).^(0:1:nterms)./factorial.(collect(big(0.):big(1.):nterms))
    pexpminusz=Polynomial(expminusz)
    logzplusone=[0; (-1).^(0:1:nterms-1)./(1:1:nterms)]
    plogzplusone=Polynomial(logzplusone)

    # From exp(z+δ) ≈ p(z), approximate δ ≈ log((exp(-z) p(z)-1)+1).
    b=Polynomial((pexpminusz * a).coeffs[1:nterms])
    c=Polynomial((plogzplusone(b-1)).coeffs[1:nterms])

    # Compute bound on |δ|.
    bnd_bwd_err=Polynomial(abs.(c.coeffs))

    # Find point where bound on relative backward error equals tolerance.
    e_bwd(z)=abs.(bnd_bwd_err(z))./abs.(z)
    h(z)=e_bwd(z) - tolerance
    theta_bwd=fzero(h, theta_init)

    return e_bwd,convert(coefftype,theta_bwd)
end
