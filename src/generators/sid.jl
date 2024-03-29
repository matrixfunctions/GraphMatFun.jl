# Boostingthecomputationofthematrixexponential

export graph_sid_exp;

"""
    (graph,cref)=graph_sid_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential using `k` matrix
multiplications following the procedure in the reference. The coefficients are
directly copied from the paper.

The evaluation is embedded in the degopt-format, see [`graph_degopt`](@ref),
using the function [`graph_sastre_yks_degopt`](@ref). Moreover, for `k<=3` it
uses the Paterson–Stockmeyer method.

**Reference**

1. J. Sastre, J. Ibáñez, E. Defez. "Boosting the computation of the matrix
   exponential". Applied Mathematics of Computation, 340:206-220, 2019.
   DOI: [10.1016/j.amc.2018.08.017](https://doi.org/10.1016/j.amc.2018.08.017)
"""
function graph_sid_exp(k; T = Float64)
    if (k <= 3)
        if k == 3
            deg = 6
        elseif k == 2
            deg = 4
        elseif k == 1
            deg = 2
        else
            deg = 1
        end
        return graph_ps_degopt(convert.(T, 1 ./ factorial.(0:deg)))

    elseif (k == 4)
        c1  = 4.018761610201036 * 1e-4
        c2  = 2.945531440279683 * 1e-3
        c3  = -8.709066576837676 * 1e-3
        c4  = 4.017568440673568 * 1e-1
        c5  = 3.230762888122312 * 1e-2
        c6  = 5.768988513026145 * 1e0
        c7  = 2.338576034271299 * 1e-2
        c8  = 2.381070373870987 * 1e-1
        c9  = 2.224209172496374 * 1e0
        c10 = -5.792361707073261 * 1e0
        c11 = -4.130276365929783 * 1e-2
        c12 = 1.040801735231354 * 1e1
        c13 = -6.331712455883370 * 1e1
        c14 = 3.484665863364574 * 1e-1
        c15 = 1
        c16 = 1

        c = [
            [[c2, c1], [0.0, 0, 0]], # (10)
            [[1.0], [0, c4, c3], [1.0], [0, 0, c5], [c6], [0.0, 0, c7]], # (11)
            [
                [0.0, 1],
                [0, c9, c8],
                [c10, 1],
                [0, c11, 0],
                [c13, c12],
                [c16, c15, c14],
            ], # (12)
        ]

        c = convert.(Vector{Vector{T}}, c)

        (graph, cref) = graph_sastre_yks_degopt(2, 2, c)

    elseif (k == 5)
        c1  = 1.161658834444880e-06
        c2  = 4.500852739573010e-06
        c3  = 5.374708803114821e-05
        c4  = 2.005403977292901e-03
        c5  = 6.974348269544424e-02
        c6  = 9.418613214806352e-01
        c7  = 2.852960512714315e-03
        c8  = -7.544837153586671e-03
        c9  = 1.829773504500424e+00
        c10 = 3.151382711608315e-02
        c11 = 1.392249143769798e-01
        c12 = -2.269101241269351e-03
        c13 = -5.394098846866402e-02
        c14 = 3.112216227982407e-01
        c15 = 9.343851261938047e+00
        c16 = 6.865706355662834e-01
        c17 = 3.233370163085380e+00
        c18 = -5.726379787260966e+00
        c19 = -1.413550099309667e-02
        c20 = -1.638413114712016e-01

        c = [
            [[c3, c2, c1], [0.0, 0, 0, 0]], # (17)
            [
                [1.0],
                [0, c6, c5, c4],
                [1.0],
                [0, 0, c8, c7],
                [c9],
                [0, 0, c11, c10],
            ], # (18)
            [
                [0.0, 1],
                [0, c14, c13, c12],
                [c15, 1],
                [0, c16, 0, 0],
                [c18, c17],
                [1, 1, c20, c19],
            ], # (19)
        ]

        c = convert.(Vector{Vector{T}}, c)

        (graph, cref) = graph_sastre_yks_degopt(2, 3, c)

    elseif (k == 6)
        # Coefficients from http://personales.upv.es/~jorsasma/software/expmpol.m
        c = [
            1.172460202011541e-08
            9.379681616092325e-08
            1.406952242413849e-06
            2.294895435403922e-05
            2.024281516007681e-03
            1.430688980356062e-02
            1.952545843107103e-01
            2.865001388641538e+00
            -1.204349003694297e-03
            2.547056607231984e-03
            2.721930992200371e-02
            2.498969092549990e+02
            2.018492049443954e-02
            1.965098904519709e-01
            1.739158441630994e+00
            8.290085751394409e+00
            2.919349464582001e-04
            1.758035313846159e-04
            1.606091400855144e-02
            3.655234395347475e-02
            2.243394407902074e-03
            -3.005000525808178e-02
            1.969779342112314e-01
            1
            1
        ]

        # Evaluation from  http://personales.upv.es/~jorsasma/software/expmpol.m
        # y0s=Ap{4}*(c(1)*Ap{4}+c(2)*Ap{3}+c(3)*Ap{2}+c(4)*A);
        # y1s=(y0s+c(5)*Ap{4}+c(6)*Ap{3}+c(7)*Ap{2}+c(8)*A)*(y0s+c(9)*Ap{4}+c(10)*Ap{3}+c(11)*Ap{2})+c(12)*y0s+c(13)*Ap{4}+c(14)*Ap{3}+c(15)*Ap{2}+c(16)*A;
        # sol=y1s*(y0s+c(17)*Ap{4}+c(18)*Ap{3}+c(19)*Ap{2}+c(20)*A)+c(21)*Ap{4}+c(22)*Ap{3}+c(23)*Ap{2}+A+eye(n);
        CC = [
            [[c[4], c[3], c[2], c[1]], [0.0, 0, 0, 0, 0]],
            [
                [1.0],
                [0, c[8], c[7], c[6], c[5]],
                [1.0],
                [0, 0, c[11], c[10], c[9]],
                [c[12]],
                [0, c[16], c[15], c[14], c[13]],
            ],
            [
                [0.0, 1],
                [0.0, 0, 0, 0, 0],
                [1.0, 0],
                [0, c[20], c[19], c[18], c[17]],
                [0.0, 0],
                [1, 1, c[23], c[22], c[21]],
            ], # (19)
        ]

        CC = convert.(Vector{Vector{T}}, CC)

        (graph, cref) = graph_sastre_yks_degopt(2, 4, CC)

    elseif (k == 7)
        # Coefficients from http://personales.upv.es/~jorsasma/software/expmpol.m
        c = [
            1.556371639324141e-11
            1.556371639324141e-10
            2.957106114715868e-09
            6.204734935438909e-08
            1.313681421698863e-06
            3.501669195497238e-05
            1.283057135586989e-03
            2.479095151834799e-02
            4.155284057336423e-01
            5.951585263506065e+00
            3.753710741641900e-05
            2.100333647757715e-04
            2.630043177655382e-03
            3.306559506631931e-02
            6.175954247606858e+01
            2.742336655922557e-03
            3.005135891320298e-02
            2.857950268422422e-01
            2.991654767354374e+00
            1.110689398085882e+01
            8.572383602707347e-06
            9.027588625491207e-05
            1.121744731945438e-03
            8.139086096860678e-03
            -2.638236222337760e-04
            6.263526066651383e-05
            4.985549176118462e-03
            7.705596948494946e-02
            5.029302610017967e-01
            1
            1
        ]

        # Evaluation from  http://personales.upv.es/~jorsasma/software/expmpol.m
        # y0s=Ap{5}*(c(1)*Ap{5}+c(2)*Ap{4}+c(3)*Ap{3}+c(4)*Ap{2}+c(5)*A);
        # y1s=(y0s+c(6)*Ap{5}+c(7)*Ap{4}+c(8)*Ap{3}+c(9)*Ap{2}+c(10)*A)*(y0s+c(11)*Ap{5}+c(12)*Ap{4}+c(13)*Ap{3}+c(14)*Ap{2})+c(15)*y0s+c(16)*Ap{5}+c(17)*Ap{4}+c(18)*Ap{3}+c(19)*Ap{2}+c(20)*A;
        # sol=y1s*(y0s+c(21)*Ap{5}+c(22)*Ap{4}+c(23)*Ap{3}+c(24)*Ap{2}+c(25)*A)+c(26)*Ap{5}+c(27)*Ap{4}+c(28)*Ap{3}+c(29)*Ap{2}+A+eye(n);
        CC = [
            [[c[5], c[4], c[3], c[2], c[1]], [0.0, 0, 0, 0, 0, 0]],
            [
                [1.0],
                [0, c[10], c[9], c[8], c[7], c[6]],
                [1.0],
                [0, 0, c[14], c[13], c[12], c[11]],
                [c[15]],
                [0, c[20], c[19], c[18], c[17], c[16]],
            ],
            [
                [0.0, 1],
                [0.0, 0, 0, 0, 0, 0],
                [1.0, 0],
                [0, c[25], c[24], c[23], c[22], c[21]],
                [0.0, 0],
                [1, 1, c[29], c[28], c[27], c[26]],
            ],
        ]

        CC = convert.(Vector{Vector{T}}, CC)

        (graph, cref) = graph_sastre_yks_degopt(2, 5, CC)

    elseif (k > 7)
        (graph, cref) = graph_sid_exp(7; T = T)
        # Square it 7-k times
        for j = 8:k
            degopt = Degopt(graph)
            square!(degopt)
            scale!(degopt, convert(T, 1 / 2))
            (graph, cref) = graph_degopt(degopt)
        end
    else
        error("Not implemented for k=$k")
    end

    return (graph, cref)
end
