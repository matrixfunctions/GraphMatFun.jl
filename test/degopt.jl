using LinearAlgebra
@testset "bbc" begin
    x1a=[3; 4.0]
    x1b=[0.1; 0.3];
    x2a=[0.1; 0.1; -0.1];
    x2b=[0.3; 0; -0.1];
    x=[(x1a,x1b), (x2a,x2b)]


    z=[1.0; -0.1; 0.1; 0.01];

    degopt1=Degopt(deepcopy(x),deepcopy(z));
    # testing gen_degopt_poly
    (g1,_)=gen_degopt_poly(degopt1)

    # Testing grow
    degopt2=grow!(deepcopy(degopt1));
    (g2,_)=gen_degopt_poly(degopt1)
    @test eval_graph(g1,0.3) ==  eval_graph(g2,0.3)


    # Test get_degopt_crefs
    (x_crefs,y_crefs)=get_degopt_crefs(2)
    set_coeffs!(g1,5.0,[x_crefs[1][2][2]]);

    x[1][2][2]=5.0;
    degopt3=Degopt(deepcopy(x),deepcopy(z));
    (g3,_)=gen_degopt_poly(degopt3)

    @test eval_graph(g1,-3.0) ==  eval_graph(g3,-3.0)


    # Test scale!

    degopt4=deepcopy(degopt1);
    scale!(degopt4,2);
    (g4,_)=gen_degopt_poly(degopt4);
    @test eval_graph(g4,0.3) == eval_graph(g1,2*0.3)


end
