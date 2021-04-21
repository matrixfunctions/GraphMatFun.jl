export opt_gauss_newton!

"""
    (iter,resnorm)=opt_gauss_newton!(graph, objfun, discr; maxit = 100, logger = 0,
                                     errtype = :abserr, stoptol = 1e-6,
                                     cref = get_all_cref(graph), γ0 = 1.0,
                                     linlsqr = :backslash, droptol = 0)

Applies the Gauss–Newton algorithm to solve the nonlinear least squares problem:
Fit the output of the `graph` to the values of `objfun`, in the points `discr`.

The variable `graph` is modified during the iterations.
The function returns the iteration `iter` where it terminated, and the corresponding residual
norm `resnorm`.

The kwargs are as follows:
`maxit` determines the maximum number of iterations.
If `logger` has a value >0, then the intermediate results are printed.
`errtype` determines the error type measured, see `adjust_for_errtype!`.
`stoptol` is the corresponding stopping tolerence.
`cref` is a `Vector{Tuple{Symbol,Int}}` that determines which coefficients of `graph`
that are considered free variables and optimized.
The stepsize can be scaled with `γ0`.
`linlsqr` and `droptol` determines how the inner linear least squares problem is
solved; see `solve_linlsqr!`.

    """
function opt_gauss_newton!(graph, objfun, discr;
                           maxit = 100,
                           logger = 0,
                           errtype = :abserr,
                           stoptol = 1e-6,
                           cref = get_all_cref(graph),
                           γ0 = 1.0,
                           linlsqr = :backslash,
                           droptol = 0
                           )

    resnorm=Inf
    iter=maxit
    vals = init_vals_eval_graph!(graph, discr, nothing)
    objfun_vals = objfun.(discr)

    for j=1:maxit
        F = eval_graph(graph, discr, vals=vals)
        res = F - objfun_vals
        Jac = eval_jac(graph, discr, cref, vals=vals)
        adjust_for_errtype!(Jac, res, objfun_vals, errtype)

        resnorm=norm(res)
        if (logger>0)
            @show j, resnorm
        end
        if (resnorm<stoptol)
            iter = j
            break
        end

        d = solve_linlsqr!(Jac, res, linlsqr, droptol)
        x = get_coeffs(graph, cref)
        x -= γ0*d
        set_coeffs!(graph, x, cref)
    end
    return (iter,resnorm)
end
