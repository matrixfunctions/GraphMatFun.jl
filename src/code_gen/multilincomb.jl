# Represents a graph which with different lincomb functionality:
# It can add more than two matrices. Only used in code generation.
struct MultiLincombCompgraph{T}
    operations::Dict{Symbol,Symbol}
    parents::Dict{Symbol,Vector{Symbol}}
    coeffs::Dict{Symbol,Vector{T}}
    outputs::Vector{Symbol}
    super_graph
end

# Workaround. This is just a wrapper
function MultiLincombCompgraph(g::Compgraph)
    T = eltype(g)
    return MultiLincombCompgraph{T}(g.operations,
                          g.parents,
                          g.coeffs,
                          g.outputs,
                          g);
end
