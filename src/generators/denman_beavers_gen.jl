export graph_denman_beavers

"""
        (graph,cref)=graph_denman_beavers(k;T=ComplexF64,cref_mode=0)

Creates the graph corresponding to the Denman-Beavers
iteration for the matrix square root, using k iterations.
The kwarg `cref_mode` specifies if references
to all X and Y (`cref_mode=0`) should be stored
or just Y (`cref_mode=-1`) or just X (`cref_mode=1`).

Reference:

* The matrix sign function and computations in systems, E. Denman,  A. Beavers, Applied Mathematics and Computation, 2(1), 1976, https://doi.org/10.1016/0096-3003(76)90020-5
    """
function graph_denman_beavers(k;T=ComplexF64,cref_mode=0)

    graph = Compgraph(T);

    cref=[];
    for j=0:k
        Xj=Symbol("X$j");
        Yj=Symbol("Y$j");

        if (j==0)
            # X0=A
            # Y0=I
            Xj=:A;
            Yj=:I;
        end

        # Yjinv=Yj^{-1}
        Yinvj=Symbol("Yinv$j");
        if (j>0)
            add_ldiv!(graph,Yinvj,Yj,:I);
        else
            Yinvj=:I; # First step Yinv0=inv(I)=I;
        end

        # Xjinv=Xj^{-1}
        add_ldiv!(graph,Symbol("Xinv$j"),Xj,:I);

        # X_{j+1}=0.5*Xj+0.5*inv(Yj)
        add_lincomb!(graph,Symbol("X"*string(j+1)),
                     0.5,Xj,0.5,Yinvj);
        if (cref_mode>=0)
            push!(cref,(Symbol("X"*string(j+1)),1))
            push!(cref,(Symbol("X"*string(j+1)),2))
        end

        # Y_{j+1}=0.5*Yj+0.5*inv(Xj)
        add_lincomb!(graph,Symbol("Y"*string(j+1)),
                     0.5,Yj,0.5,Symbol("Xinv$j"));
        if (cref_mode<=0)
            push!(cref,(Symbol("Y"*string(j+1)),1))
            push!(cref,(Symbol("Y"*string(j+1)),2))
        end

        add_output!(graph,Symbol("X"*string(j+1)));
    end
    # Not needed to get the output
    del_node!(graph,Symbol("Y"*string(k+1)))
    setdiff!(cref,[(Symbol("Y"*string(k+1)),1),(Symbol("Y"*string(k+1)),2)])

    del_node!(graph,Symbol("Xinv"*string(k)));
#
    return (graph,cref);
end
