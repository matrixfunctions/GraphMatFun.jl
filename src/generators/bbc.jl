export gen_general_poly_recursion


"""
    (graph,crefs)=gen_general_poly_recursion(k;compress_keys=true,T=ComplexF64)
    (graph,crefs)=gen_general_poly_recursion(x,z;compress_keys=true)

Corresponds to the general polynomial recursion

    B1=A
    B2=(x *I+x *A)(x *I+x *A)
    B3=(x *I+x *A+x*B2)(x *I+x *A+x *B2)
     ..

and

    Out=z*I+z*A1+z*B2+z*B3+z*B4..

The `x`-values are given in the argument `x`, which is
a `Vector{Tuple{Vector,Vector}}`, containing the elements
of each sum. The `z`-vector contains the elements
to form the output. If `compress_keys=true`, the references
to `z[3],z[4],...` are not returned.
If the parameter `k` is supplied instead of the coefficients,
all coeffs will be set to one.

Reference: The general recursion is mentioned in equation (9) in this paper:

* Philipp Bader, Sergio Blanes, and Fernando Casas. Computing the matrix exponential with an optimized Taylor polynomial approximation. Mathematics, 7(12), 2019.
"""
function gen_general_poly_recursion(x,z;compress_keys=true)
    (graph,crefs)=gen_general_poly_recursion_B(x);

    k=size(x,1);

    # Add the polynomial in the Bk coeffs

    z_nodes=[:I; :A];
    for s=2:k+1
        push!(z_nodes,Symbol("B$s"));
    end
    key=Symbol("T2k$(k+3)");
    crefs_new=add_sum!(graph,key,
                       z,z_nodes,Symbol("T2k"));


    if (compress_keys)
        # Only take the I and A
        append!(crefs,crefs_new[1:2]);
    else
        append!(crefs,crefs_new);
    end


    # Set the output
    empty!(graph.outputs);
    push!(graph.outputs,key);

    return (graph,crefs);

end
function gen_general_poly_recursion(k;T=ComplexF64,compress_keys=true)

    x=Vector{Tuple{Vector{T},Vector{T}}}(undef,k);
    for i=1:k
        x[i]=(ones(T,i+1),ones(T,i+1));
    end

    z=ones(T,k+2)

    gen_general_poly_recursion(x,z;compress_keys=compress_keys)
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


        base0="B$s";
        base_a="Ba$s";
        base_b="Bb$s";

        # First poly
        c_a=x[s-1][1]
        crefs_a=add_sum!(graph,Symbol(base_a),c_a,useful_syms,
                         Symbol("$(base_a)_"));

        append!(cref,crefs_a);

        # Second poly
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