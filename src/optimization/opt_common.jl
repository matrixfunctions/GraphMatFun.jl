export adjust_for_errtype!, solve_linlsqr

"""
    adjust_for_errtype!(A, b, objfun_vals, errtype)

Adjusts the matrix `A` and vector `b` to `errtype`.
`A` is a matrix, e.g., the Jacobian or system matrix in a least squares problem.
`b` is a vector, e.g., the residuals or the right-hand side in a lieast squares problem.
`objfun_vals` is a vector of objective function values, and
`errtype` can be `:abserr` or `:relerr`, i.e., absolute error or relative error.
     """
function adjust_for_errtype!(A, b, objfun_vals, errtype)
    # Assumes a Jacobian and residual-vector that is computed in absolute error.
    # Adjusts Jacobian and residual-vector from absolute, to relative error, i.e.,
    # objective function 1/2 sum_i (Z(x_i)-f(x_i))^2 / f(x_i)^2 and
    if (errtype==:abserr)
        # Do nothing
    elseif (errtype==:relerr)
        objfun_vals[objfun_vals.==0] .= eps()*100 # Hack to avoid division by zero
        D = Diagonal(1 ./ objfun_vals)
        A[:] =D*A
        b[:] =D*b
    else
        error("Unknown errtype '", errtype, "'.")
    end
    return (A, b)
end


"""
    d = solve_linlsqr(A, b, linlsqr, droptol)

Solves the linear least squares problem

    Ad=b.

`linlsqr` determines how the linear least squares problem is solved.
It can be `:backslash`, `:nrmeq`, or `:svd`, and for the latter
singular values below `droptol` are disregarded.
     """
function solve_linlsqr(A, b, linlsqr, droptol)
    if (linlsqr == :backslash)
        d = A\b
    elseif (linlsqr == :nrmeq)
        d = (A'*A)\(A'*b)
    elseif (linlsqr == :svd)
        if (eltype(A) == BigFloat || eltype(A) == Complex{BigFloat})
            # You must use include "using GenericSVD"
            Sfact=svd!(A; full=false, alg=nothing)
        else
            Sfact=svd(A)
        end
        d=Sfact.S
        # Use pseudoinverse if droptol>0
        II = (d/d[1]) .< droptol
        dinv = 1 ./ d
        dinv[II] .= 0
        # No explicit construction, only multiplication
        # JJ0=Sfact.U*Diagonal(d)*Sfact.Vt
        d = Sfact.V*((Diagonal(dinv))*(Sfact.U'*b));
    end
    return d
end
