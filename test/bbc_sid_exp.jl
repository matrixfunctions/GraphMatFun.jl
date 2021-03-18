using LinearAlgebra
@testset "bbc sid exp" begin

    α=big(1.1);

    x=1.1;
    p2=4; # Order of approx for 2 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(2,T=BigFloat);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)


    @test -log(abs(err2/err1))/log(α) > p2+1



    x=0.8;
    p3=8; # Order of approx for 3 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(3);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)

    @test -log(abs(err2/err1))/log(α) > p3+1


    x=2.0
    p4bbc=12; # Order of approx for 4 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(4);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)

    @test -log(abs(err2/err1))/log(α) > p4bbc+1


    p4sid=15; # Order of approx for 4 matrix multiplies
    (graph,cref)=gen_sid_exp(4);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p4sid+1


    x=2.3
    p5bbc=18 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p5bbc+1


    x=3.2;
    p5sid=21 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_sid_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p5sid+1




    x=3.5;
    p6sid=24 # Order of approx for 6 matrix multiplies
    (graph,cref)=gen_sid_exp(6);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p6sid+1


    x=5.2;
    p7sid=30 # Order of approx for 7 matrix multiplies
    (graph,cref)=gen_sid_exp(7);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p7sid+1



end
