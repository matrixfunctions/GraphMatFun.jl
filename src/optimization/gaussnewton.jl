export opt_gaussnewton!

"""
    (iter,resnorm)=opt_gaussnewton!(graph, objfun, discr; maxit = 100, logger = 0,
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
`errtype` determines the error type measured, see `adjust_for_errtype()`, and
`stoptol` is the corresponding stopping tolerence.
`cref` is a `Vector{Tuple{Symbol,Int}}` that determines which coefficients of `graph`
that are considered free variables and optimized.
The stepsize can be scaled with `γ0`.
`linlsqr` determines how the inner linear least squares problem is solved.
It can be `:backslash`, `:nrmeq`, or `:svd`, and for the latter
singular values below `droptol` are disregarded.


    """
function opt_gaussnewton!(graph, objfun, discr;
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

    for j=1:maxit
        F = eval_graph(graph, discr, vals=vals)
        objfun_vals = objfun.(discr)
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

        d = solve_linlsqr(Jac, res, linlsqr, droptol)
        x = get_coeffs(graph, cref)
        x += γ0*d
        set_coeffs!(graph, x, cref)
    end
    return (iter,resnorm)
end


#Internal to sovle the linear least squares problem in GN
function solve_linlsqr(Jac, res, linlsqr, droptol)
    if (linlsqr == :backslash)
        d = -Jac\res
    elseif (linlsqr == :nrmeq)
        d = -(Jac'*Jac)\(Jac'*res)
    elseif (linlsqr == :svd)
        if (eltype(Jac) == BigFloat || eltype(Jac) == Complex{BigFloat})
            # You must use include "using GenericSVD"
            Sfact=svd!(Jac; full=false, alg=nothing)
        else
            Sfact=svd(Jac)
        end
        d=Sfact.S
        # Use pseudoinverse if droptol>0
        II = (d/d[1]) .< droptol
        dinv = 1 ./ d
        dinv[II] .= 0
        # No explicit construction, only multiplication
        # JJ0=Sfact.U*Diagonal(d)*Sfact.Vt
        d=-Sfact.V*((Diagonal(dinv))*(Sfact.U'*res));
    end
end
