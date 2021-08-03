export graph_bbcs_cheb_exp

"""
    (graph,cref)=graph_bbcs_cheb_exp(k;T=Complex{BigFloat})

Computes a polynomial evaluation approximating the exponential for
skew-Hermitian matrices using `k` matrix multiplications following the procedure
in the reference.

For `k` > 5 it resorts to scaling-and-squaring of the `k` = 5 graph.
The graph is in the degopt format, see [`graph_degopt`](@ref).

Reference:

[^BBCS21]: P. Bader, S. Blanes, F. Casas, M. Seydaoglu. "An efficient algorithm
    to compute the exponential of skew-Hermitian matrices for the time integration
    of the Schrödinger equation".
    [arXiv:2103.10132 [math.NA]](https://arxiv.org/abs/2103.10132), 2021.
"""
function graph_bbcs_cheb_exp(k; T = Complex{BigFloat})
    CBF = Complex{BigFloat}
    if (k <= 5) # Basic implementation from paper
        if (k == 1)
            x1a = vcat(zero(CBF), one(CBF))
            x1b = vcat(zero(CBF), one(CBF))
            xv = [(x1a, x1b)]
            y = [
                BigFloat("0.9999999999999999999998"),
                BigFloat("-0.9999999999761950000001")im,
                BigFloat("-0.4999999999920650000000"),
            ]
        elseif (k == 2)
            x1a = vcat(zero(CBF), one(CBF))
            x1b = vcat(zero(CBF), one(CBF))
            x2a = vcat(zeros(CBF, 2), one(CBF))
            x2b = vcat(
                zero(CBF),
                BigFloat("0.16666657785001893215")im,
                BigFloat("0.04166664890333648869"),
            )
            xv = [(x1a, x1b), (x2a, x2b)]
            y = [
                BigFloat("0.99999999999999999997"),
                BigFloat("-0.99999999999981067844")im,
                BigFloat("-0.49999999999994320353"),
                one(CBF),
            ]
        elseif (k == 3) # Equation (22)
            x1a = vcat(zero(CBF), one(CBF))
            x1b = vcat(zero(CBF), one(CBF))
            x2a = vcat(zeros(CBF, 2), one(CBF))
            x2b = vcat(
                zero(CBF),
                BigFloat("431") / BigFloat("4000"),
                BigFloat("-0.02693906873598870733")im,
            )
            x3a = vcat(
                zeros(CBF, 2),
                BigFloat("0.66321004441662438593")im,
                one(CBF),
            )
            x3b = vcat(
                BigFloat("0.54960853911436015786")im,
                BigFloat("0.16200952846773660904"),
                BigFloat("-0.01417981805211804396")im,
                BigFloat("-0.03415953916892111403"),
            )
            xv = [(x1a, x1b), (x2a, x2b), (x3a, x3b)]
            y = [
                BigFloat("0.99999999999999999928"),
                BigFloat("-0.99999999999999233987")im,
                BigFloat("-0.13549409636220703066"),
                zero(CBF),
                one(CBF),
            ]
        elseif (k == 4) # Equation (23)
            x1a = vcat(zero(CBF), one(CBF))
            x1b = vcat(zero(CBF), one(CBF))
            x2a = vcat(zeros(CBF, 2), one(CBF))
            x2b = vcat(zero(CBF), one(CBF), zero(CBF))
            # B_4^2
            x3a = vcat(
                zero(CBF),
                BigFloat("0.13340427306445612526")im,
                BigFloat("0.02022602029818310774"),
                BigFloat("-0.00674638241111650999")im,
            )
            x3b = x3a
            # (B_2 + A_6)*A_6 = (B_2 + B_3 + B_4^2)*(B_3 + B_4^2)
            b2 = vcat(
                zero(CBF),
                BigFloat("1.41183797496250375498")im,
                zero(CBF),
                BigFloat("-0.00866935318616372016")im,
            )
            b3 = vcat(
                BigFloat("2.69584306915332564689"),
                BigFloat("-1.35910926168869260391")im,
                BigFloat("-0.09896214548845831754"),
                BigFloat("0.01596479463299466666")im,
            )
            x4a = vcat((b2 + b3), one(CBF))
            x4b = vcat(b3, one(CBF))
            xv = [(x1a, x1b), (x2a, x2b), (x3a, x3b), (x4a, x4b)]
            # B_1 + (B_2 + A_6)*A_6 = B_1 + previous
            b1 = vcat(
                BigFloat("-6.26756985350202252845"),
                BigFloat("2.52179694712098096140")im,
                BigFloat("0.05786296656487001838"),
                BigFloat("-0.07766686408071870344")im,
            )
            y = vcat(b1, zero(CBF), one(CBF))
        elseif (k == 5) # Equation (25)
            x1a = vcat(zero(CBF), one(CBF))
            x1b = vcat(zero(CBF), one(CBF))
            x2a = vcat(zeros(CBF, 2), one(CBF))
            x2b = vcat(zero(CBF), one(CBF), zero(CBF))
            x3a = vcat(zeros(CBF, 3), one(CBF))
            x3b = vcat(zeros(CBF, 3), one(CBF))
            # B_1 * B_5
            x4a = vcat(
                zero(CBF),
                BigFloat("3") / BigFloat("25"),
                BigFloat("-0.00877476096879703859")im,
                BigFloat("-0.00097848453523780954"),
                zero(CBF),
            )
            x4b = vcat(
                zeros(CBF, 2),
                BigFloat("-0.123953695858283131480")im,
                BigFloat("-0.011202694841085592373"),
                BigFloat("-0.000012367240538259896")im,
            )

            # (B_3 + A_9)*A_9 = (B_3 + B_4 + previous)*(B_4 + previos)
            b3 = vcat(
                BigFloat("-2.58175430371188142440"),
                BigFloat("-1.73033278310812419209")im,
                BigFloat("-0.07673476833423340755"),
                BigFloat("-0.00261502969893897079")im,
                BigFloat("-0.00003400011993049304"),
            )
            b4 = vcat(
                BigFloat("2.92377758396553673559"),
                BigFloat("1.44513300347488268510")im,
                BigFloat("0.12408183566550450221"),
                BigFloat("-0.01957157093642723948")im,
                BigFloat("0.00002425253007433925"),
            )
            x5a = vcat((b3 + b4), one(CBF))
            x5b = vcat(b4, one(CBF))
            xv = [(x1a, x1b), (x2a, x2b), (x3a, x3b), (x4a, x4b), (x5a, x5b)]
            # B_2 + (B_3 + A_9)*A_9 = B_2 + previous
            b2 = vcat(
                zero(CBF),
                BigFloat("-0.66040840760771318751")im,
                BigFloat("-1.09302278471564897987"),
                BigFloat("0.25377155817710873323")im,
                BigFloat("0.00054374267434731225"),
            )
            y = vcat(b2, zero(CBF), one(CBF))
        end
        # Force convert to type
        xv = map(
            i -> (convert.(T, xv[i][1]), convert.(T, xv[i][2])),
            1:size(xv, 1),
        )
        y = convert.(T, y)
        return graph_degopt(xv, y)
    else  # Scalign-and-squaring phase
        s = k - 5
        degopt = Degopt(graph_bbcs_cheb_exp(5; T = T)[1])
        scale!(degopt, BigFloat("2.0")^(-big(s)))
        for i = 1:s
            square!(degopt)
        end
        return graph_degopt(degopt)
    end
end
