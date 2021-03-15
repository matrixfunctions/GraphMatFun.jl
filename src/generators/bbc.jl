export gen_general_poly_recursion


"""
    (graph,crefs)=gen_bbc(k;compress_keys=true,T=ComplexF64)

Corresponds to the general polynomial recursion

  B1=A
  B2=(x1*I+x2*A)(x3*I+x4*A)
  B3=(x4*I+x5*A+x6*B2)(x7*I+x8*A+x8*B2)
   ..

and

  Out=xk*I+x(k+1)*A1+B2+B3+B4..

If `compress_keys=true`, then cref contains `k+1`
keys corresponding to `x1`,... `x(k+1)`.


Reference: The general recursion is described in this paper:

* Philipp Bader, Sergio Blanes, and Fernando Casas. Computing the matrix exponential with an optimized Taylor polynomial approximation. Mathematics, 7(12), 2019.
"""
function gen_general_poly_recursion(k;T=ComplexF64,compress_keys=true)

    (graph,cref)=gen_general_poly_recursion_B(k,T);

    # Add the polynomial in the Bk coeffs
    key=:T2k1
    add_lincomb!(graph,key,1.0,:I,1.0,:A);
    push!(cref,(key,1))
    push!(cref,(key,2))

    for s=2:k+1
        prevkey=key;
        key=Symbol("T2k$s");
        add_lincomb!(graph,key,1.0,prevkey,1.0,Symbol("B$s"));
        if (!compress_keys)
            # Corresponds to a free variable in Bk so could
            # be removed
            push!(cref,(key,2));
        end
    end

    empty!(graph.outputs);
    push!(graph.outputs,key);

    return (graph,cref);

end

#
function gen_general_poly_recursion_B(k,T)
    x=Vector{Tuple{Vector{T},Vector{T}}}(undef,k);
    for i=1:k
        x[i]=(ones(T,i+1),ones(T,i+1));
    end
    return gen_general_poly_recursion_B(x);
end



# Normally x::Vector{Tuple{Vector{Number},Vector{Number}}}
# Containing the coefficients in the recursion
function gen_general_poly_recursion_B(x)


    k=size(x,1);
    T=eltype(eltype(eltype(x)));

    graph=Compgraph(T);
    useful_syms=[:I,:A];

    cref=[];
    for s=2:k+1

        #@show useful_syms;
        base0="B$s";
        base_a="Ba$s";
        base_b="Bb$s";

        # First poly
        c_a=x[s-1][1]
        crefs_a=add_sum!(graph,Symbol(base_a),c_a,useful_syms,
                         Symbol("$(base_a)_"));

        append!(cref,crefs_a);

        c_b=x[s-1][2]
        crefs_b=add_sum!(graph,Symbol(base_b),c_b,useful_syms,
                         Symbol("$(base_b)_"));
        append!(cref,crefs_b);


        # Multiply them together
        add_mult!(graph,Symbol(base0),Symbol(base_a),Symbol(base_b));
        push!(useful_syms,Symbol(base0))
        empty!(graph.outputs);
        push!(graph.outputs,Symbol(base0));
    end

    return (graph,cref);
end



function gen_general_poly_recursion_B_old(k,T)

    graph=Compgraph(T);
    useful_syms=[:I,:A];

    cref=[];
    for s=2:k+1

        #@show useful_syms;
        base0="B$s";
        base1="Ba$s";
        base2="Bb$s";
        # First poly
        base=base1
        if (size(useful_syms,1)==2)
            key=Symbol("$(base)");
        else
            key=Symbol("$(base)_2");
        end

        push!(cref,(key,1))
        push!(cref,(key,2))
        add_lincomb!(graph,key,1.0,useful_syms[1],1.0,useful_syms[2])

        for (j,termkey)=enumerate(useful_syms[3:end]);
            prevkey=key;
            if (termkey != useful_syms[end])
                key=Symbol("$(base)_$(j+2)");
            else
                key=Symbol(base);
            end
            push!(cref,(key,2))
            add_lincomb!(graph,key,1.0,prevkey,1.0,termkey)
        end


        # Second poly
        base=base2;
        if (size(useful_syms,1)==2)
            key=Symbol("$(base)");
        else
            key=Symbol("$(base)_2");
        end
        push!(cref,(key,1))
        push!(cref,(key,2))
        add_lincomb!(graph,key,1.0,useful_syms[1],1.0,useful_syms[2])

        for (j,termkey)=enumerate(useful_syms[3:end]);
            prevkey=key;
            if (termkey != useful_syms[end])
                key=Symbol("$(base)_$(j+2)");
            else
                key=Symbol(base);
            end

            push!(cref,(key,2))
            add_lincomb!(graph,key,1.0,prevkey,1.0,termkey)
        end

        # Multiply them together
        add_mult!(graph,Symbol(base0),Symbol(base1),Symbol(base2));
        push!(useful_syms,Symbol(base0))
        empty!(graph.outputs);
        push!(graph.outputs,Symbol(base0));
    end

    return (graph,cref);
end
