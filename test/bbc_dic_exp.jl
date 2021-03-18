using LinearAlgebra
@testset "bbc dic exp" begin

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


    p4dic=15; # Order of approx for 4 matrix multiplies
    (graph,cref)=gen_dic_exp(4);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p4dic+1


    x=2.3
    p5bbc=18 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_bbc_basic_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p5bbc+1


    x=3.2;
    p5dic=21 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_dic_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p5dic+1

end
