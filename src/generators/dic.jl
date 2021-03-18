# Boostingthecomputationofthematrixexponential

export gen_dic_exp;

"""
    (graph,cref)=gen_dic_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference. The coefficients are directly copied from the paper.

Reference:

* Boosting the computation of the matrix exponential, Sastre, Ibanez, Defez, Appl. Math. Computation, 340, 2019, 206-220.

    """
function gen_dic_exp(k;T=Float64)

    if (k==4)

        c1=convert(T,4.018761610201036*1e-4)
        c2=convert(T,2.945531440279683*1e-3)
        c3=convert(T,-8.709066576837676*1e-3)
        c4=convert(T,4.017568440673568*1e-1)
        c5=convert(T,3.230762888122312*1e-2)
        c6=convert(T,5.768988513026145*1e0)
        c7=convert(T,2.338576034271299*1e-2)
        c8=convert(T,2.381070373870987*1e-1)
        c9=convert(T,2.224209172496374*1e0)
        c10=convert(T,-5.792361707073261*1e0)
        c11=convert(T,-4.130276365929783*1e-2)
        c12=convert(T,1.040801735231354*1e1)
        c13=convert(T,-6.331712455883370*1e1)
        c14=convert(T,3.484665863364574*1e-1)
        c15=convert(T,1)
        c16=convert(T,1)


        # B2=A^2
        v1a=[0;1.0];
        v1b=[0;1.0];

        # B3=y02=A^2*(c1*A^2+c2*A)
        v2a=[0;0;1];
        v2b=[0;c2;c1];

        # B4=(c4*A+c3*A^2+y02)*(c5*A^2+y02)   # First term in y12
        #       Note: y12=B4+c6*y02+c7*A^2
        v3a=[0;c4;c3;1.0];
        v3b=[0;0;c5;1.0];


        # B5=(c9*A+c8*A^2+y12)*(c11*A+c10*y02+y12)    # first term in y22
        #   =(c9*A+(c8+c7)*A^2+c6*y02+B4)*(c11*A+c7*A^2+(c10+c6)*y02+B4)
        v4a=[0;c9;c8+c7;c6;1.0];
        v4b=[0;c11;c7;(c10+c6);1.0];


        # add the additional terms to B5

        y=[c16;c15;(c14+c12*c7);(c13+c12*c6);c12;1.0];

        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_general_poly_recursion(xv,y);


    elseif (k==5)
        c=[1.161658834444880e-06
           4.500852739573010e-06
           5.374708803114821e-05
           2.005403977292901e-03
           6.974348269544424e-02
           9.418613214806352e-01
           2.852960512714315e-03
           -7.544837153586671e-03
           1.829773504500424e+00
           3.151382711608315e-02
           1.392249143769798e-01
           -2.269101241269351e-03
           -5.394098846866402e-02
           3.112216227982407e-01
           9.343851261938047e+00
           6.865706355662834e-01
           3.233370163085380e+00
           -5.726379787260966e+00
           -1.413550099309667e-02
           -1.638413114712016e-01
           1
           1];

        c=convert.(T,c);


        # B2= A^2
        v1a=[0;1.0];
        v1b=[0;1.0];

        # B3= A^3
        v2a=[0;1.0;0];
        v2b=[0;0;1.0];

        # B4=y03
        v3a=[0;0;0;1]; # A^3
        v3b=[0;c[3];c[2];c[1]];

        # B5=first term in y13
        #   = (c6*A+c[5]*A^2+c[4]*A^3+y03)*(c8*A^2+c7*A^3+y03)
        # Note: y13=B5+c11*A^2+c10*A^3+c9*y03
        v4a=[0;c[6];c[5];c[4];1];
        v4b=[0;0;c[8];c[7];1];

        y13=[0;0;c[11];c[10];c[9];1];


        # B6=first term in y23
        #   =(c14*A+c13*A^2+      c12*A^3+          y13)*(c16*A+c15*y03+y13)
        #   =(c14*A+(c13+c11)*A^2+(c12+c10)*A^3+c9*y03+B5)*
        #               (c16*A+c11*A^2+c10*A^3+(c15+c9)*y03+B5)

        #v5a=[0;c[14];(c[13]+c[11]);(c[12]+c[10]);c[9];1];
        #v5b=[0;c[16];c[11];c[10];(c[15]+c[9]);1];
        v5a=[0;c[14];c[13];c[12];0;0]+y13;
        v5b=[0;c[16];0;0;c[15];0]+y13;


        y=[1;1;c[20];c[19];c[18];0;1];
        y[1:6] += c[17]*y13;


        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b); (v5a,v5b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_general_poly_recursion(xv,y);


    else

        error("Not implemented for k=$k");
    end


    return (graph,cref);
end
