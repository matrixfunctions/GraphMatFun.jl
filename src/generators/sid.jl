# Boostingthecomputationofthematrixexponential

export gen_sid_exp;

"""
    (graph,cref)=gen_sid_exp(k;T=Float64)

Computes a polynomial evaluation approximating the exponential
using `k` matrix multiplications following the procedure
in the reference. The coefficients are directly copied from the paper.

Reference:

* Boosting the computation of the matrix exponential, Sastre, Ibanez, Defez, Appl. Math. Computation, 340, 2019, 206-220.

    """
function gen_sid_exp(k;T=Float64)

    if (k==4)

        c=[
            4.018761610201036*1e-4
            2.945531440279683*1e-3
            -8.709066576837676*1e-3
            4.017568440673568*1e-1
            3.230762888122312*1e-2
            5.768988513026145*1e0
            2.338576034271299*1e-2
            2.381070373870987*1e-1
            2.224209172496374*1e0
            -5.792361707073261*1e0
            -4.130276365929783*1e-2
            1.040801735231354*1e1
            -6.331712455883370*1e1
            3.484665863364574*1e-1
            1
            1]

        c=convert.(T,c);

        # B2=A^2
        v1a=[0;one(T)];
        v1b=[0;one(T)];

        # B3=y02=A^2*(c1*A^2+c2*A)
        v2a=[0;0;1];
        v2b=[0;c[2];c[1]];

        # B4=(c4*A+c3*A^2+y02)*(c5*A^2+y02)   # First term in y12
        #       Note: y12=B4+c6*y02+c7*A^2
        v3a=[0;c[4];c[3];1.0];
        v3b=[0;0;c[5];1.0];


        # B5=(c9*A+c8*A^2+y12)*(c11*A+c10*y02+y12)    # first term in y22
        #   =(c9*A+(c8+c7)*A^2+c6*y02+B4)*(c11*A+c7*A^2+(c10+c6)*y02+B4)
        v4a=[0;c[9];c[8]+c[7];c[6];1.0];
        v4b=[0;c[11];c[7];(c[10]+c[6]);1.0];


        # add the additional terms to B5

        y=[c[16];c[15];(c[14]+c[12]*c[7]);(c[13]+c[12]*c[6]);c[12];1.0];

        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_degopt_poly(xv,y);


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

        (graph,cref)=gen_degopt_poly(xv,y);


    elseif (k==6)
        # Coefficients from http://personales.upv.es/~jorsasma/software/expmpol.m
        c=[1.172460202011541e-08
           9.379681616092325e-08
           1.406952242413849e-06
           2.294895435403922e-05
           2.024281516007681e-03
           1.430688980356062e-02
           1.952545843107103e-01
           2.865001388641538e+00
           -1.204349003694297e-03
           2.547056607231984e-03
           2.721930992200371e-02
           2.498969092549990e+02
           2.018492049443954e-02
           1.965098904519709e-01
           1.739158441630994e+00
           8.290085751394409e+00
           2.919349464582001e-04
           1.758035313846159e-04
           1.606091400855144e-02
           3.655234395347475e-02
           2.243394407902074e-03
           -3.005000525808178e-02
           1.969779342112314e-01
           1
           1];
        c=convert.(T,c);


        # B2= A^2
        v1a=[0;1.0];
        v1b=[0;1.0];

        # B3= A^3
        v2a=[0;1.0;0];
        v2b=[0;0;1.0];

        # B4= A^4
        v3a=[0;0;1.0;0];
        v3b=[0;0;1.0;0];


        # B5=y04=A^4*(c4*A+c3*A^2+c2*A^3+c1*A^4)
        v4a=[0;0;0;0;1.0];
        v4b=[0;c[4];c[3];c[2];c[1]];

        # B6=first term in y14
        #   =(c8*A+c7*A^2+c6*A^3+c5*A^4+y04)*(c11*A^2+c10*A^3+c9*A^4+y04)
        v5a=[0;c[8];c[7];c[6];c[5];1];
        v5b=[0;0;c[11];c[10];c[9];1];

        # The additional terms in y14
        y14=[0;c[16];c[15];c[14];c[13];c[12];1];




        # B7=first term in T24
        #   = y14*(c20*A+c19*A^2+c18*A^3+c17*A^4+y04)
        v6a=y14;
        v6b=[0;c[20];c[19];c[18];c[17];1;0];

        # Output
        y=[1;1;c[23];c[22];c[21];0;0;1];

        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b); (v5a,v5b); (v6a,v6b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_degopt_poly(xv,y);


    elseif (k==7)
        # Coefficients from http://personales.upv.es/~jorsasma/software/expmpol.m
        c=[1.556371639324141e-11
           1.556371639324141e-10
           2.957106114715868e-09
           6.204734935438909e-08
           1.313681421698863e-06
           3.501669195497238e-05
           1.283057135586989e-03
           2.479095151834799e-02
           4.155284057336423e-01
           5.951585263506065e+00
           3.753710741641900e-05
           2.100333647757715e-04
           2.630043177655382e-03
           3.306559506631931e-02
           6.175954247606858e+01
           2.742336655922557e-03
           3.005135891320298e-02
           2.857950268422422e-01
           2.991654767354374e+00
           1.110689398085882e+01
           8.572383602707347e-06
           9.027588625491207e-05
           1.121744731945438e-03
           8.139086096860678e-03
           -2.638236222337760e-04
           6.263526066651383e-05
           4.985549176118462e-03
           7.705596948494946e-02
           5.029302610017967e-01
           1
           1];

        c=convert.(T,c);

        # B2= A^2
        v1a=[0;1.0];
        v1b=[0;1.0];

        # B3= A^3
        v2a=[0;1.0;0];
        v2b=[0;0;1.0];

        # B4= A^4
        v3a=[0;0;1.0;0];
        v3b=[0;0;1.0;0];

        # B5= A^5
        v4a=[0;1.0;0.0;0;0];
        v4b=[0;0;  0.0;0;1.0];


        # B6=y05=
        #   A^5*(c5*A+c4*A^2+c3*A^3+c2*A^4+c1*A^5)
        v5a=[0;0.0; 0.0; 0;  0;    1];
        v5b=[0;c[5];c[4];c[3];c[2];c[1]];

        # B7=first term in y15
        #   =(c10*A+c9*A^2+c8*A^3+c7*A^4+c6*A^5+y05)*
        #         (c14*A^2+c13*A^3+c12*A^4+c11*A^5+y05)
        v6a=[0;c[10];c[9] ;c[8] ;c[7]; c[6];1]
        v6b=[0;0    ;c[14];c[13];c[12];c[11];1]

        y15=[0;c[20];c[19];c[18];c[17];c[16];c[15];1];

        # B8=first term in T30
        #   =y15*(c[25]*A+c[24]*A^2+c[23]*A^3+c[22]*A^4+c[21]*A^5+y05)
        v7a=y15;
        v7b=[0;c[25];c[24];c[23];c[22];c[21];1;0];

        # Output
        y=[1;1;c[29];c[28];c[27];c[26];0;0;1];


        xv=[(v1a,v1b); (v2a,v2b); (v3a,v3b); (v4a,v4b); (v5a,v5b); (v6a,v6b); (v7a,v7b)];

        # Force convert to type
        xv=map(i-> (convert.(T,xv[i][1]),convert.(T,xv[i][2])),1:size(xv,1))
        y = convert.(T,y);

        (graph,cref)=gen_degopt_poly(xv,y);
    elseif (k>7)
        (graph,cref)=gen_sid_exp(7;T=T)
        # Square it 7-k times
        for j=8:k

            degopt=Degopt(graph);
            square!(degopt);
            scale!(degopt,convert(T,1/2));
            (graph,cref)=gen_degopt_poly(degopt);
        end
    else
        error("Not implemented for k=$k");
    end


    return (graph,cref);
end
