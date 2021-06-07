# Represents a graph which with different lincomb functionality:
# It can add more than two matrices. Only used in code generation.
struct MultiLincombCompgraph{T}
    operations::Dict{Symbol,Symbol}
    parents::Dict{Symbol,NTuple{<:Any,Symbol}}
    coeffs::Dict{Symbol,NTuple{<:Any,T}}
    outputs::Vector{Symbol}
end
function MultiLincombCompgraph(g::Compgraph)
    T=eltype(g);
    extgraph=MultiLincombCompgraph(Dict{Symbol,Symbol}(),Dict{Symbol,NTuple{<:Any,Symbol}}(),Dict{Symbol,NTuple{<:Any,T}}(),Vector{Symbol}())
    Z=extract_sums(g);
    for k in keys(g.operations)
        if (g.operations[k] == :mult)
            add_mult!(extgraph,k,g.parents[k][1],g.parents[k][2]);
        end
        if (g.operations[k] == :ldiv)
            add_ldiv!(extgraph,k,g.parents[k][1],g.parents[k][2]);
        end
    end
    for s in Z
        coeff_list=s[1];
        symbol_list=s[2];
        key=s[3][end];
        p=size(coeff_list,1);
        extgraph.operations[key]=:lincomb;
        extgraph.coeffs[key]=NTuple{p,T}(coeff_list);
        extgraph.parents[key]=NTuple{p,Symbol}(symbol_list);
    end
    for k in g.outputs;
        push!(extgraph.outputs,k);
    end
    return extgraph

end
