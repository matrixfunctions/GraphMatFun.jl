using LinearAlgebra
@testset "bbc exp" begin

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
    p4=12; # Order of approx for 4 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(4);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)


    @test -log(abs(err2/err1))/log(α) > p4+1


    x=2.3
    p5=18 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p5+1

end
