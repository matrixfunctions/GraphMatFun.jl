export gen_bbcs_cheb_exp

"""
    (graph,cref)=gen_bbcs_cheb_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference.

Reference:

*  An efficient algorithm to compute the exponential of skew-Hermitian matrices for the time integration of the Schrödinger equation, P. Bader, S. Blanes, F. Casas, M. Seydaoglu, arXiv:2103.10132
    """
function gen_bbcs_cheb_exp(k;T=ComplexF64)
    if (k==1)
        x1a = vcat(zero(T), one(T))
        x1b = vcat(zero(T), one(T))
        xv = [ (x1a, x1b) ]
        y = [0.9999999999999999999998, -0.9999999999761950000001im, -0.4999999999920650000000]
    elseif (k==2)
        x1a = vcat(zero(T), one(T))
        x1b = vcat(zero(T), one(T))
        x2a = vcat(zeros(T,2), one(T))
        x2b = vcat(zero(T), 0.16666657785001893215im, 0.04166664890333648869)
        xv = [ (x1a, x1b),
               (x2a, x2b) ]
        y = [0.99999999999999999997, -0.99999999999981067844im, -0.49999999999994320353, one(T)]
    elseif (k==3)
        x1a = vcat(zero(T), one(T))
        x1b = vcat(zero(T), one(T))
        x2a = vcat(zeros(T,2), one(T))
        x2b = vcat(zero(T), 431/4000, -0.02693906873598870733im)
        x3a = vcat(zeros(T,2), 0.66321004441662438593im, one(T))
        x3b = vcat(0.54960853911436015786im, 0.16200952846773660904, -0.01417981805211804396im, -0.03415953916892111403)
        xv = [ (x1a, x1b),
               (x2a, x2b),
               (x3a, x3b) ]
        y = [0.99999999999999999928, -0.99999999999999233987im, -0.13549409636220703066, zero(T), one(T)]
    elseif (k==4)
        x1a = vcat(zero(T), one(T))
        x1b = vcat(zero(T), one(T))
        x2a = vcat(zeros(T,2), one(T))
        x2b = vcat(zero(T), one(T), zero(T))
        # B_4^2
        x3a = vcat(zero(T), 0.13340427306445612526im, 0.02022602029818310774, -0.00674638241111650999im)
        x3b = x3a
        # (B_2 + A_6)*A_6 = (B_2 + B_3 + B_4^2)*(B_3 + B_4^2)
        b2 = [zero(T), 1.41183797496250375498im, zero(T), -0.00866935318616372016im]
        b3 = [2.69584306915332564689, -1.35910926168869260391im, -0.09896214548845831754, 0.01596479463299466666im]
        x4a = vcat((b2+b3), one(T))
        x4b = vcat(b3, one(T))
        xv = [ (x1a, x1b),
               (x2a, x2b),
               (x3a, x3b),
               (x4a, x4b) ]
        # B_1 + (B_2 + A_6)*A_6 = B_1 + previous
        b1 = vcat(-6.26756985350202252845, 2.52179694712098096140im, 0.05786296656487001838, -0.07766686408071870344im)
        y = vcat(b1, zero(T), one(T))
    else
        error("Not implemented k=$k")
    end

    # Force convert to type
    xv = map( i -> (convert.(T,xv[i][1]),convert.(T,xv[i][2])), 1:size(xv,1) )
    y = convert.(T,y)
    return gen_degopt_poly(xv,y);
end
