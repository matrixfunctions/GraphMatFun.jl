using LinearAlgebra
@testset "sastre" begin

    x=2.0
    p1=12; # Order of approx for 4 matrix multiplies
    (graph,cref)=gen_sastre_basic_exp(4);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    α=big(1.1);
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)


    @test -log(abs(err2/err1))/log(α) > p1+1


    x=2.3
    p2=16 # Order of approx for 5 matrix multiplies
    (graph,cref)=gen_sastre_basic_exp(5);
    err1=eval_graph(big(graph),big(x))-exp(big(x))
    err2=eval_graph(big(graph),big(x)/α)-exp(big(x)/α)
    @test -log(abs(err2/err1))/log(α) > p1+1

end
