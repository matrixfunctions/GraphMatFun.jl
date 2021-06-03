using LinearAlgebra
@testset "bbc" begin


    x1a=[3; 4.0]
    x1b=[0.1; 0.3];
    x2a=[0.1; 0.1; -0.1];
    x2b=[0.3; 0; -0.1];
    x=[(x1a,x1b), (x2a,x2b)]


    z=[1.0; -0.1; 0.1; 0.01];

    (graph,crefs)=gen_degopt_poly(x,z,compress_keys=false)

    x0=3.1415;
    ev1=eval_graph(graph,x0,comporder=get_topo_order_degopt(2));

    # Manual evaluation
    β1=(x1a[1]+x1a[2]*x0)*(x1b[1]+x1b[2]*x0)
    β2=(x2a[1]+x2a[2]*x0+x2a[3]*β1)*(x2b[1]+x2b[2]*x0+x2b[3]*β1)
    ev2=z[1]+z[2]*x0+z[3]*β1+z[4]*β2;

    @test ev1≈ev2


    # Also check that the crefs are correct

    # Modify the graph: set almost all coeffs to one
    T=Float64;
    t=ones(T,size(crefs,1));
    t[2]=2.0;
    set_coeffs!(graph,t,crefs);

    # Create a new graph object graph2: with coeffs almost all zero
    x1a=[1; 2.0]
    x1b=[1; 1.0];
    x2a=[1.0; 1.0; 1.0];
    x2b=[1.0; 1.0; 1.0];

    z=ones(T,4);
    x=[(x1a,x1b), (x2a,x2b)]
    (graph2,crefs2)=gen_degopt_poly(x,z)



    x1=0.3;

    @test eval_graph(graph,x1,comporder=get_topo_order_degopt(2))==eval_graph(graph2,x1);


    for k = 1:7
        x=Vector{Tuple{Vector{T},Vector{T}}}(undef,k)
        for i=1:k
            x[i]=(collect(1.0:i+1),collect(-1.0:-1:-(i+1)))
        end
        z=collect(range(4,5,length=k+2))
        (graph,cref) = gen_degopt_poly(x,z)

        (xx,zz) = get_degopt_crefs(k)
        for i=1:k
                for N = 1:2
                    @test get_coeffs(graph,xx[i][N]) == x[i][N]
                end
        end
        @test get_coeffs(graph,zz) == z
    end


end
