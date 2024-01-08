@testset "bigraph" begin
    x=[0
       -1.7157450924607922
       0.8507824740023241
       0.8142424866036422
       -0.5044212203692009
       -0.44309270492554415
       -0.9261617361310625
       -0.9908932075112057
       0.540893328611887
       -0.6896982802624144]
    # Create something with both many terms and a single term
    (graph0,_)=graph_ps_degopt(x);
    compress_graph_zero_coeff!(graph0)
    graph0.coeffs[:Bb2][1]=-0.9;
    rename_node!(graph0,:Ba5,:Bb42) # Improve test coverage
    x=[3.0 4; 1.1 0.1];
    (graph1,_)=graph_bigraph(graph0)
    @test eval_graph(graph0,x) == eval_graph(graph1,x)
end
