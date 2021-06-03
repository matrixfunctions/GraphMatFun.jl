export gen_bbcs_cheb_exp

"""
    (graph,cref)=gen_bbcs_cheb_exp(k;T=Complex{BigFloat})

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference.

Reference:

*  An efficient algorithm to compute the exponential of skew-Hermitian matrices for the time integration of the Schrödinger equation, P. Bader, S. Blanes, F. Casas, M. Seydaoglu, arXiv:2103.10132
    """
function gen_bbcs_cheb_exp(k;T=Complex{BigFloat})
    CBF = Complex{BigFloat}
    if (k==1)
        x1a = vcat(zero(CBF), one(CBF))
        x1b = vcat(zero(CBF), one(CBF))
        xv = [ (x1a, x1b) ]
        y = [BigFloat("0.9999999999999999999998"), BigFloat("-0.9999999999761950000001")im, BigFloat("-0.4999999999920650000000")]
    elseif (k==2)
        x1a = vcat(zero(CBF), one(CBF))
        x1b = vcat(zero(CBF), one(CBF))
        x2a = vcat(zeros(CBF,2), one(CBF))
        x2b = vcat(zero(CBF), BigFloat("0.16666657785001893215")im, BigFloat("0.04166664890333648869"))
        xv = [ (x1a, x1b),
               (x2a, x2b) ]
        y = [BigFloat("0.99999999999999999997"), BigFloat("-0.99999999999981067844")im, BigFloat("-0.49999999999994320353"), one(CBF)]
    elseif (k==3)
        x1a = vcat(zero(CBF), one(CBF))
        x1b = vcat(zero(CBF), one(CBF))
        x2a = vcat(zeros(CBF,2), one(CBF))
        x2b = vcat(zero(CBF), BigFloat("431")/BigFloat("4000"), BigFloat("-0.02693906873598870733")im)
        x3a = vcat(zeros(CBF,2), BigFloat("0.66321004441662438593")im, one(CBF))
        x3b = vcat(BigFloat("0.54960853911436015786")im, BigFloat("0.16200952846773660904"), BigFloat("-0.01417981805211804396")im, BigFloat("-0.03415953916892111403"))
        xv = [ (x1a, x1b),
               (x2a, x2b),
               (x3a, x3b) ]
        y = [BigFloat("0.99999999999999999928"), BigFloat("-0.99999999999999233987")im, BigFloat("-0.13549409636220703066"), zero(CBF), one(CBF)]
    elseif (k==4)
        x1a = vcat(zero(CBF), one(CBF))
        x1b = vcat(zero(CBF), one(CBF))
        x2a = vcat(zeros(CBF,2), one(CBF))
        x2b = vcat(zero(CBF), one(CBF), zero(CBF))
        # B_4^2
        x3a = vcat(zero(CBF), BigFloat("0.13340427306445612526")im, BigFloat("0.02022602029818310774"), BigFloat("-0.00674638241111650999")im)
        x3b = x3a
        # (B_2 + A_6)*A_6 = (B_2 + B_3 + B_4^2)*(B_3 + B_4^2)
        b2 = vcat(zero(CBF), BigFloat("1.41183797496250375498")im, zero(CBF), BigFloat("-0.00866935318616372016")im)
        b3 = vcat(BigFloat("2.69584306915332564689"), BigFloat("-1.35910926168869260391")im, BigFloat("-0.09896214548845831754"), BigFloat("0.01596479463299466666")im)
        x4a = vcat((b2+b3), one(CBF))
        x4b = vcat(b3, one(CBF))
        xv = [ (x1a, x1b),
               (x2a, x2b),
               (x3a, x3b),
               (x4a, x4b) ]
        # B_1 + (B_2 + A_6)*A_6 = B_1 + previous
        b1 = vcat(BigFloat("-6.26756985350202252845"), BigFloat("2.52179694712098096140")im, BigFloat("0.05786296656487001838"), BigFloat("-0.07766686408071870344")im)
        y = vcat(b1, zero(CBF), one(CBF))
    else
        error("Not implemented k=$k")
    end

    # Force convert to type
    xv = map( i -> (convert.(T,xv[i][1]),convert.(T,xv[i][2])), 1:size(xv,1) )
    y = convert.(T,y)
    return gen_degopt_poly(xv,y);
end
