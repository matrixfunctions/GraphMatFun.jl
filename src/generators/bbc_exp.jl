export graph_bbc_exp

"""
    (graph,cref)=graph_bbc_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential using `k` matrix
multiplications following the procedure in the reference. The coefficients are
directly copied from the paper.

The evaluation is embedded in the `degopt`-format, see [`graph_degopt`](@ref),
and for `k`<3, the evaluation is using the Paterson–Stockmeyer method.

Reference:

[1] P. Bader, S. Blanes, and F. Casas. "Computing the matrix exponential with an
    optimized Taylor polynomial approximation". Mathematics, 7(12), 2019. DOI:
    [10.3390/math7121174](https://doi.org/10.3390/math7121174)
"""
function graph_bbc_exp(k; T = Float64)
    if (k < 3)
        if k == 2
            deg = 4
        elseif k == 1
            deg = 2
        else
            deg = 1
        end
        return graph_ps_degopt(convert.(T, 1 ./ factorial.(0:deg)))

    elseif (k == 3)
        # Equation 13
        value177 = convert(T, 177) # Lessen roundoff errors
        y0 = 1
        y1 = 1
        y2 = (857 - 58 * sqrt(value177)) / 630
        x3 = 2 / 3
        x1 = x3 * (1 + sqrt(value177)) / 88
        x2 = x3 * (1 + sqrt(value177)) / 352
        x4 = (-271 + 29 * sqrt(value177)) / (315 * x3)
        x5 = 11 * (-1 + sqrt(value177)) / (1260 * x3)
        x6 = 11 * (-9 + sqrt(value177)) / (5040 * x3)
        x7 = (89 - sqrt(value177)) / (5040 * x3^2)

        # A2=A^2
        v1a = [0; 1.0]
        v1b = [0; 1.0]

        # A4= A2 * (x1*A+x2*A2)
        v2a = [0; 0; 1.0]
        v2b = [0; x1; x2]

        # A8=(x3*A2+A4)*(x4*I+x5*A+x6*A2+x7*A4)
        v3a = [0; 0; x3; 1.0]
        v3b = [x4; x5; x6; x7]

        y = [y0; y1; y2; 0; 1.0]

        xv = [(v1a, v1b); (v2a, v2b); (v3a, v3b)]
        (graph, cref) = graph_degopt(xv, y)

    elseif (k == 4)
        a01 = -0.01860232051462055322
        a02 = 4.60000000000000000000
        a03 = 0.21169311829980944294
        a04 = 0
        a11 = -0.00500702322573317730
        a12 = 0.99287510353848683614
        a13 = 0.15822438471572672537
        a14 = -0.13181061013830184015
        a21 = -0.57342012296052226390
        a22 = -0.13244556105279963884
        a23 = 0.16563516943672741501
        a24 = -0.02027855540589259079
        a31 = -0.13339969394389205970
        a32 = 0.00172990000000000000
        a33 = 0.01078627793157924250
        a34 = -0.00675951846863086359

        v1a = [0; 1.0]
        v1b = [0; 1.0]
        v2a = [0; 0; 1.0]
        v2b = [0; 1; 0.0]

        # B4^2
        v3a = [a04; a14; a24; a34]
        v3b = [a04; a14; a24; a34]

        # (B2+A6)*A6
        v4b = [a03; a13; a23; a33; 1] # A6
        v4a = v4b + [a02; a12; a22; a32; 0] # B2+A6

        y = [a01; a11; a21; a31; 0; 1]

        #xv=[v1,v2,v3,v4];
        xv = [(v1a, v1b); (v2a, v2b); (v3a, v3b); (v4a, v4b)]
        (graph, cref) = graph_degopt(xv, y)

    elseif (k == 5)
        a01 = 0
        a11 = -0.10036558103014462001
        a21 = -0.00802924648241156960
        a31 = -0.00089213849804572995

        b01 = 0
        b11 = 0.39784974949964507614
        b21 = 1.36783778460411719922
        b31 = 0.49828962252538267755
        b61 = -0.00063789819459472330

        b02 = -10.9676396052962062593
        b12 = 1.68015813878906197182
        b22 = 0.05717798464788655127
        b32 = -0.00698210122488052084
        b62 = 0.00003349750170860705

        b03 = -0.09043168323908105619
        b13 = -0.06764045190713819075
        b23 = 0.06759613017704596460
        b33 = 0.02955525704293155274
        b63 = -0.00001391802575160607

        b04 = 0
        b14 = 0
        b24 = -0.09233646193671185927
        b34 = -0.01693649390020817171
        b64 = -0.00001400867981820361

        v1a = [0; 1.0]
        v1b = [0; 1.0]
        v2a = [0; 1.0; 0.0]
        v2b = [0; 0; 1.0]
        v3a = [0; 0; 0; 1.0]
        v3b = [0; 0; 0; 1.0]

        # B1B5
        v4a = [a01; a11; a21; a31; 0]
        v4b = [b04; b14; b24; b34; b64]
        # (B3+A9) * A9
        v5b = [b03; b13; b23; b33; b63; 1]  # A9
        v5a = v5b + [b02; b12; b22; b32; b62; 0]  # A9+B3
        v5 = vcat(v5a, v5b)

        # T18: B2+(B3+A9)*A9
        y = [b01; b11; b21; b31; b61; 0; 1]

        xv = [(v1a, v1b); (v2a, v2b); (v3a, v3b); (v4a, v4b); (v5a, v5b)]
        (graph, cref) = graph_degopt(xv, y)

    elseif (k == 6)
        # Coefficients unknown
        error("k=6 for BBC the coefficients are not documented.")
    else
        error("Incorrect k")
    end
    return (graph, cref)
end
