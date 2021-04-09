using LinearAlgebra
@testset "rational" begin

    A = [-3 1 ; 5 6.6];

    # Simple calls with coefficients
    den_coeff = [3; -1; 2.0; 0.1; 1.2; 0.9; -0.3; 1.6];
    num_coeff = [1.4, 2.1, -4, -0.1, 0.2, -1.2, 2.1, 0.4, 1.4];
    den = den_coeff[1]*I + den_coeff[2]*A + den_coeff[3]*A^2 + den_coeff[4]*A^3 + den_coeff[5]*A^4 +
          den_coeff[6]*A^5 + den_coeff[7]*A^6 + den_coeff[8]*A^7
    num = num_coeff[1]*I + num_coeff[2]*A + num_coeff[3]*A^2 + num_coeff[4]*A^3 + num_coeff[5]*A^4 +
          num_coeff[6]*A^5 + num_coeff[7]*A^6 + num_coeff[8]*A^7 + num_coeff[9]*A^8
    P = num\den

    (graph, cref) = gen_rational(den_coeff, num_coeff, gen_monomial)
    @test eval_graph(graph,A) ≈ P

    (graph, cref) = gen_rational(den_coeff, num_coeff, gen_ps)
    @test eval_graph(graph,A) ≈ P


    # Generate more complicated graphs before
    x1a=[1.2; 4.4]
    x1b=[0.3; -0.3]
    x2a=[2.1; 0.1; -0.1]
    x2b=[0.3; 0.2; 0.1]
    x=[(x1a,x1b), (x2a,x2b)]
    z=[1.0; 0.1; -0.1; 0.01]
    den_coeffs=[(x1a,x1b), (x2a,x2b)]
    β1=(x1a[1]*I+x1a[2]*A)*(x1b[1]*I+x1b[2]*A)
    β2=(x2a[1]*I+x2a[2]*A+x2a[3]*β1)*(x2b[1]*I+x2b[2]*A+x2b[3]*β1)
    den=z[1]*I+z[2]*A+z[3]*β1+z[4]*β2;
    (den_graph,den_cref)=gen_degopt_poly(x,z,compress_keys=false)

    x1a=[2.1; 1.2]
    x1b=[1.3; 0.3]
    x2a=[0.1; -1; 1]
    x2b=[0; 2; -0.3]
    x=[(x1a,x1b), (x2a,x2b)]
    z=[1.0; -1; 0.1; -0.01]
    num_coeffs=[(x1a,x1b), (x2a,x2b)]
    β1=(x1a[1]*I+x1a[2]*A)*(x1b[1]*I+x1b[2]*A)
    β2=(x2a[1]*I+x2a[2]*A+x2a[3]*β1)*(x2b[1]*I+x2b[2]*A+x2b[3]*β1)
    num=z[1]*I+z[2]*A+z[3]*β1+z[4]*β2;
    (num_graph,num_cref)=gen_degopt_poly(x,z,compress_keys=false)

    P = num\den
    (graph, cref) = gen_rational(den_graph, num_graph, den_cref=den_cref, num_cref=num_cref)
    @test eval_graph(graph,A) ≈ P

end
