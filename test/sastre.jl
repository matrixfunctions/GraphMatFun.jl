using Polynomials

@testset "Sastre" begin


    ## Test approximation of exponential and number of multiplications
    for (k,Ak,m) = [(0,0.0001,:ps_degopt), (1,0.001,:ps_degopt), (2,0.1,:ps_degopt), (3,1.1,:y1s), (4,2.9,:y1s), (6,5.1,:h2m), (8,16.5,:z1ps)] # Where do I find these limits?
        A = randn(100,100)/20 * Ak
        for i = 1:2
            if i == 1
                (graph,cref) = graph_sastre_exp(k,m)
            else
                (graph,cref) = graph_sastre_exp(k,:auto)
            end
            @test eval_graph(graph,A) ≈ exp(A)
            @test sum(values(graph.operations) .== :mult) == k
        end
    end



    ## Test the degree-9 solver graph_sastre_poly(b) by comparison to tabulated values
    # 1: Table 4 against equations (16)-(32) for exponential
    (graph1,cref1) = graph_sastre_exp(3,:y1s)
    (graph2,cref2) =  graph_sastre_poly(1 ./factorial.(0:8))
    @test sum(values(graph2.operations) .== :mult) == 3
    @test all(cref1 .== cref2)
    for c = cref1
        @test graph1.coeffs[c[1]][c[2]] ≈ graph2.coeffs[c[1]][c[2]]
    end

    # 2: Table 11 against Table 10 and equations (16)-(32)
    (graph_exp,_) = graph_sastre_exp(6,:h2m)
    A = randn(100,100)/20 * 5.3
    # Naive implementation from coefficients in Table 10
    b8  = 2.186201576339059e-7
    bp8 = 2.186201576339059e-7
    b7  = 9.839057366529322e-7
    bp7 = 2.514016785489562e-6
    b6  = 1.058964584814256e-5
    bp6 = 3.056479369585950e-5
    b5  = 1.554700173279057e-4
    bp5 = 3.197607034851565e-4
    b4  = 2.256892506343887e-3
    bp4 = 2.585006547542889e-3
    b3  = 2.358987357109499e-2
    bp3 = 1.619043970183846e-2
    b2  = 1.673139636901279e-1
    bp2 = 8.092036376147299e-2
    b1  = 7.723603212944010e-1
    bp1 = 3.229486011362677e-1
    b0  = 3.096467971936040
    bp0 = 0.0
    β0 = 1
    (graph,_) =  graph_sastre_poly([b0,b1,b2,b3,b4,b5,b6,b7,b8])
    (graphp,_) =  graph_sastre_poly([bp0,bp1,bp2,bp3,bp4,bp5,bp6,bp7,bp8])
    v1 = eval_graph(graph,A)
    vp1 = eval_graph(graphp,A)

    @test v1*vp1 + β0*I ≈ eval_graph(graph_exp,A)
    @test v1*vp1 + β0*I ≈ exp(A)



    ## Test conversion of (62)-(65) to Degopt
    # 1: By using numbers from (57)-(59)
    c16 =  4.018761610201036e-4
    c8  =  2.116367017255747e0
    c15 =  2.945531440279683e-3
    c7  = -5.792361707073261e0
    c14 =  8.712167566050691e-2
    c6  = -1.491449188999246e-1
    c13 =  4.017568440673568e-1
    c5  =  1.040801735231354e1
    c12 = -6.352311335612147e-2
    c4  = -6.331712455883370e1
    c11 =  2.684264296504340e-1
    c3  =  3.484665863364574e-1
    c10 =  1.857143141426026e1
    c2  = -1.224230230553340e-1
    c9  =  2.381070373870987e-1
    c1  =  1

    c = [
    [[c15, c16], [0.0, 0, 0]], # (57)
    [[1.0], [0, c13, c14], [1.0], [c11, 0, c12], [c10], [0.0, 0, 0]], # (58)
    [[0.0, 1], [0, c8, c9], [c7, 1], [0, c6, 0], [c4, c5], [c1, c2, c3]], # (59)
    ]

    (graph,cref)=graph_sastre_yks_degopt(2,2,c)
    A = randn(100,100)/20 * 3.1
    @test eval_graph(graph,A) ≈ exp(A)

    x = Polynomial("x")
    CC = coeffs(eval_graph(Compgraph(Any,graph),x))
    t = (CC[17]-1/factorial(big(16)))*factorial(big(16)) -(-0.454) # (60)
    @test (abs(t) < 1e-3) && (abs(t) > 1e-4) # (60)
    @test all(abs.((CC .* factorial.(big.(0:16)))[1:16] .- 1) .<= 2e-15)

    # 2: By constructing an onw polynomial of the form y_33
    c1  = 1.1e-2;    c2  = 1.2e-2;    c3  = 1.3e-2;    c4  = 1.4e-2
    c5  = 1.5e-2;    c6  = 1.6e-2;    c7  = 1.7e-2;    c8  = 2.05e-2
    c9  = 2.1e-2;    c10 = 2.15e-2;   c11 = 2.2e-2;    c12 = 2.25e-2
    c13 = 2.3e-2;    c14 = 2.35e-2;   c15 = 2.4e-2;    c16 = 2.45e-2
    c17 = 2.5e-2;    c18 = 2.55e-2;   c19 = 2.6e-2;    c20 = 2.65e-2
    c21 = 2.7e-2;    c22 = 2.75e-2;   c23 = 3.05e-2;   c24 = 3.1e-2
    c25 = 3.15e-2;   c26 = 3.2e-2;    c27 = 3.25e-2;   c28 = 3.3e-2
    c29 = 3.35e-2;   c30 = 3.4e-2;    c31 = 3.45e-2;   c32 = 3.5e-2
    c33 = 3.55e-2;   c34 = 3.6e-2;    c35 = 3.65e-2;   c36 = 3.7e-2
    c37 = 3.75e-2;   c38 = 3.8e-2;    c39 = 3.85e-2;   c40 = 3.9e-2
    c41 = 4.05e-2;   c42 = 4.1e-2;    c43 = 4.15e-2;   c44 = 4.2e-2
    c45 = 4.25e-2;   c46 = 4.3e-2;    c47 = 4.35e-2;   c48 = 4.4e-2
    c49 = 4.45e-2;   c50 = 4.5e-2;    c51 = 4.55e-2;   c52 = 4.6e-2
    c53 = 4.65e-2;   c54 = 4.7e-2;    c55 = 4.8e-2;    c56 = 4.85e-2
    c57 = 4.9e-2;    c58 = 4.95e-2;   c59 = 5e-2;     c60 = 5.05e-2
    c61 = 5.10e-2

    c = [
    [[c1, c2, c3], [c4, c5, c6, c7]],
    [[c8], [c9, c10, c11, c12], [c13], [c14, c15, c16, c17], [c18], [c19, c20, c21, c22]],
    [[c23, c24], [c25, c26, c27, c28], [c29, c30], [c31, c32, c33, c34], [c35, c36], [c37, c38, c39, c40]],
    [[c41, c42, c43], [c44, c45, c46, c47], [c48, c49, c50], [c51, c52, c53, c54], [c55, c56, c57], [c58, c59, c60, c61]],
    ]

    II = Matrix(1.0I,100,100)
    A2 = A*A
    A3 = A2*A
    y03 = A3*(c1*A+c2*A2+c3*A3) + c4*II+c5*A+c6*A2+c7*A3
    (graph,cref)=graph_sastre_yks_degopt(0,3,c[1:1])
    @test eval_graph(graph,A) ≈ y03

    y13 = (c8*y03+c9*II+c10*A+c11*A2+c12*A3)*(c13*y03+c14*II+c15*A+c16*A2+c17*A3) +
          c18*y03+c19*II+c20*A+c21*A2+c22*A3
    (graph,cref)=graph_sastre_yks_degopt(1,3,c[1:2])
    @test eval_graph(graph,A) ≈ y13

    y23 = (c23*y03+c24*y13+c25*II+c26*A+c27*A2+c28*A3)*(c29*y03+c30*y13+c31*II+c32*A+c33*A2+c34*A3) +
          c35*y03+c36*y13+c37*II+c38*A+c39*A2+c40*A3
    (graph,cref)=graph_sastre_yks_degopt(2,3,c[1:3])
    @test eval_graph(graph,A) ≈ y23

    y33 = (c41*y03+c42*y13+c43*y23+c44*II+c45*A+c46*A2+c47*A3)*(c48*y03+c49*y13+c50*y23+c51*II+c52*A+c53*A2+c54*A3) +
          c55*y03+c56*y13+c57*y23+c58*II+c59*A+c60*A2+c61*A3
    (graph,cref)=graph_sastre_yks_degopt(3,3,c)
    @test eval_graph(graph,A) ≈ y33

end
