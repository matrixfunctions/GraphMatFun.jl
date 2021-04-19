export gen_degopt_poly, get_degopt_coeffs, get_topo_order_degopt


"""
    (graph,crefs)=gen_degopt_poly(k;compress_keys=true,T=ComplexF64,input=:A)
    (graph,crefs)=gen_degopt_poly(x,z;compress_keys=true,input=:A)

Corresponds to the (for a fixed numer of multiplications) degree-optimal polynomial

    B1=A
    B2=(x *I+x *A)(x *I+x *A)
    B3=(x *I+x *A+x*B2)(x *I+x *A+x *B2)
     ..

and

    Out=z*I+z*A1+z*B2+z*B3+z*B4...

The `x`-values are given in the argument `x`, which is
a `Vector{Tuple{Vector,Vector}}`, containing the elements
of each sum. The `z`-vector contains the elements
to form the output, and `input` determines the name of the matrix A above.
If `compress_keys=true`, the references
to `z[3],z[4],...` are not returned.
If the parameter `k` is supplied instead of the coefficients,
all coeffs will be set to one.

Reference: The general recursion is mentioned in equation (9) in this paper:

* Philipp Bader, Sergio Blanes, and Fernando Casas. Computing the matrix exponential with an optimized Taylor polynomial approximation. Mathematics, 7(12), 2019.
"""
function gen_degopt_poly(x,z;compress_keys=true,input=:A)
    T = promote_type(eltype(eltype(eltype(x))), eltype(z))
    (graph,crefs)=gen_degopt_poly_B(x,T,input=input);

    k=size(x,1);

    # Add the polynomial in the Bk coeffs

    z_nodes=[:I; input];
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
function gen_degopt_poly(k;T=ComplexF64,compress_keys=true,input=:A)

    x=Vector{Tuple{Vector{T},Vector{T}}}(undef,k);
    for i=1:k
        x[i]=(ones(T,i+1),ones(T,i+1));
    end

    z=ones(T,k+2)

    gen_degopt_poly(x,z;compress_keys=compress_keys,input=:A)
end


# Normally x::Vector{Tuple{Vector{Number},Vector{Number}}}
# Containing the coefficients in the recursion
function gen_degopt_poly_B(x,T;input=:A)


    k=size(x,1);

    graph=Compgraph(T);
    useful_syms=[:I,input];

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


"""
    (x,z)=get_degopt_coeffs(k;compress_keys=false)

Retruns crefs related to `gen_degopt_poly`. Specifically
`x` is a `Vector{Tuple{Vector{Tuple{Symbol,Int}},Vector{Tuple{Symbol,Int}}}}`
such that `x[2][1]` corresponds to the coefficients of the left hand side of the
multiplication

    B2=(α_2_1 *I + α_2_2 *A)(β_2_1 *I + β_2_2 *A)

i.e., the crefs corresponding to `[α_2_1, α_2_2]`. See `gen_degopt_poly`. Hence,
`get_coeffs(graph,x[2][1])` returns the corresponding numerical values of the coefficients.

"""
function get_degopt_coeffs(k;compress_keys=false)

    TT = Tuple{Symbol,Int}
    x=Vector{Tuple{Vector{TT},Vector{TT}}}(undef,k)
    for i = 1:k
        x[i] = ( Vector{TT}(undef,i+1), Vector{TT}(undef,i+1) )
    end

    for s=2:k+1

        for (b,base) = enumerate(["Ba$s", "Bb$s"])
            if s == 2
                x[s-1][b][1] = (Symbol(base), 1)
            end
            for i = 2:s-1
                if i == 2
                    x[s-1][b][1] = (Symbol("$(base)_2"), 1)
                end
                x[s-1][b][i] = (Symbol("$(base)_$(i)"), 2)
            end
            x[s-1][b][s] = (Symbol(base), 2)
        end
    end

    zlength = compress_keys ? 2 : k+2
    z = Vector{TT}(undef,zlength)
    z[1] = (:T2k2, 1)
    z[2] = (:T2k2, 2)
    for i = 3:zlength-1
        z[i] = ( Symbol("T2k$(i)"), 2 )
    end
    if !compress_keys
        z[k+2] = ( Symbol("T2k$(k+3)"), 2 )
    end

    return (x,z)
end

"""
    order=get_topo_order_degopt(k)

A special implementation of `get_topo_order` for degree optimal polynomials
generated with `gen_degopt_poly`. The natural order of computation is to compute
row by row. See also `get_degopt_coeffs`.

"""
function get_topo_order_degopt(k)
    (x,z) = get_degopt_coeffs(k,compress_keys=false)
    computation_order=Vector{Symbol}(undef, 0)
    for i = 1:k
        for n = 1:2
            for j = 2:(i+1)
                    push!(computation_order, x[i][n][j][1])
            end
        end
        push!(computation_order, Symbol("B$(i+1)"))
    end
    for i = 2:(2+k)
        push!(computation_order, z[i][1])
    end
    return computation_order
end
