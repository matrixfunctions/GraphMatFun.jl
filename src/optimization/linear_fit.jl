export opt_linear_fit!

"""
    opt_linear_fit!(graph, objfun, discr, linear_cref;
                    errtype = :abserr,
                    linlsqr = :backslash,
                    droptol = 0)

Linear fitting of a `graph` of the form

    c_1 g_1(x) + c_2 g_2(x) + … + c_n g_n(x)

to the values of `objfun`, in the points `discr`. Reference to the coefficients
`c_1,…,c_n` should be given `linear_cref`.

The variable `graph` is modified during the iterations and the function has no
return value.

See `opt_gaussnewton!` for a description the kwarg `errtype`,
and `solve_linlsqr` for the kwargs `linlsqr`, and `droptol`.

    """
function opt_linear_fit!(graph, objfun, discr, linear_cref;
                         errtype = :abserr,
                         linlsqr = :backslash,
                         droptol = 0)

    objfun_vals = objfun.(discr)

    vals = init_vals_eval_graph!(graph, discr, nothing)
    eval_graph(graph, discr, vals=vals)

    n = length(linear_cref)
    T = eltype(valtype(vals))
    A = ones(T, size(discr,1), n)
    for k=1:n
        parent = graph.parents[linear_cref[k][1]][linear_cref[k][2]];
        A[:,k] = vals[parent]
    end

    adjust_for_errtype!(A, objfun_vals, objfun_vals, errtype)

    c = solve_linlsqr(A, objfun_vals, linlsqr, droptol)
    set_coeffs!(graph, c, linear_cref)

    return nothing
end




# function get_coeffs_linear_fit(discr,f,n0=size(discr,1);structure=:none,errtype=:abserr)
#     T=eltype(discr);
#     n=n0;
#     if (structure==:even)
#         n=Int(round(n0/2));
#         if (mod(n0,2)==1)
#             n=n+1;
#         end
#         multby=discr.^2
#     else
#         multby=discr;
#     end
#     A=ones(T,size(discr,1),n)
#
#
#     for k=2:n
#         A[:,k]=multby.*A[:,k-1];
#     end
#
#     bb=f;
#     if (!(f isa AbstractVector))
#         bb=f.(discr);
#     end
#     if (errtype==:relerr)
#         A=Diagonal(bb)\A;
#         bb=ones(eltype(bb),size(bb,1));
#     end
#
#     c=A\bb;
#
#     if (structure==:even)
#         c_even=c;
#         c=zeros(T,n0);
#         c[1:2:end]=c_even;
#     end
#
#
#     return c;
# end
