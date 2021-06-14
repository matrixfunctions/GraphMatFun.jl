using Polynomials

@testset "Sastre" begin


    ## Test approximation of exponential and number of multiplications
    for (k,Ak,m) = [(3,1.1,:y1s), (4,2.9,:y1s), (6,5.1,:h2m), (8,17.3,:z1ps)] # Where do I find these limits?
        A = randn(100,100)/20 * Ak
        (graph,cref) = gen_sastre_basic_exp(k,m)
        @test eval_graph(graph,A) ≈ exp(A)
        @test sum(values(graph.operations) .== :mult) == k
    end



    ## Test the degree-9 solver gen_sastre_basic(b) by comparison to tabulated values
    # 1: Table 4 against equations (16)-(32) for exponential
    (graph1,cref1) = gen_sastre_basic_exp(3,:y1s)
    (graph2,cref2) =  gen_sastre_basic(1 ./factorial.(0:8))
    @test sum(values(graph2.operations) .== :mult) == 3
    @test all(cref1 .== cref2)
    for c = cref1
        @test graph1.coeffs[c[1]][c[2]] ≈ graph2.coeffs[c[1]][c[2]]
    end

    # 2: Table 11 against Table 10 and equations (16)-(32)
    (graph_exp,_) = gen_sastre_basic_exp(6,:h2m)
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
    (graph,_) =  gen_sastre_basic([b0,b1,b2,b3,b4,b5,b6,b7,b8])
    (graphp,_) =  gen_sastre_basic([bp0,bp1,bp2,bp3,bp4,bp5,bp6,bp7,bp8])
    v1 = eval_graph(graph,A)
    vp1 = eval_graph(graphp,A)

    @test v1*vp1 + β0*I ≈ eval_graph(graph_exp,A)
    @test v1*vp1 + β0*I ≈ exp(A)



    ## Test conversion of (62)-(65) by using numbers from (57)-(59)
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

    (graph,cref)=sastre_yks_to_degopt(2,2,c)
    A = randn(100,100)/20 * 3.1
    @test eval_graph(graph,A) ≈ exp(A)

    x = Polynomial("x")
    CC = coeffs(eval_graph(Compgraph(Any,graph),x))
    t = (CC[17]-1/factorial(16))*factorial(16) -(-0.454) # (60)
    @test (abs(t) < 1e-3) && (abs(t) > 1e-4) # (60)
    @test all(abs.((CC .* factorial.(0:16))[1:16] .- 1) .<= 2e-15)
end
