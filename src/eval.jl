export eval_graph, eval_jac, eval_runerr

# init_vals_eval_graph!()
# init_relerrs_eval_runerr!()

## Evaluating the graph

# eval_graph()
# Modifies the Dict vals if passed as a kwarg
# Returns the output specificed in graph.outputs[output]
"""
    result=eval_graph(
        graph,
        x;
        vals = nothing,
        input = :A,
        output = size(graph.outputs, 1),
        comporder = nothing,
    )

Evaluates a graph in the value `x` which is typically a scalar value or a
matrix. If `x` is a Vector, the values will be evaluated elementwise.

The `comporder` is a `Vector` of nodes specifying in which order the graph
should be computed. By default [`get_topo_order`](@ref) is used.

The `output` is an `Int` specifying which node should be considered as output.
The output node is `graph.outputs[output]`.

The `vals` is used to inspect contents other than the output inside the graph.
Typically `vals` is a `Dict`. It will be modified to contain the computed nodes
of the graph. If we wish to inspect node `X3`, we can initiate an empty dict as
input:

```julia-repl
julia > vals = Dict{Symbol,Any}();
julia > eval_graph(graph, A, vals = vals);
julia > vals[:X3]
```
"""
function eval_graph(
    graph,
    x;
    vals = nothing,
    input = :A,
    output = size(graph.outputs, 1),
    comporder = nothing,
)
    vals = init_vals_eval_graph!(graph, x, vals, input)
    if comporder == nothing
        comporder = get_topo_order(graph, input = input)[1]
    end
    for node in comporder
        parentval1 = vals[graph.parents[node][1]]
        parentval2 = vals[graph.parents[node][2]]
        carry_out!(graph, vals, (parentval1, parentval2), node)
    end
    return vals[graph.outputs[output]]
end

# carry_out! -- Perform the operations. vectors interpreted elementwise
function carry_out!(graph, vals, parentvals, node)

    # Generic type: Treat as "scalar"
    op = graph.operations[node]
    if (op == :mult)
        if (parentsvals[1] isa Vector)
            vals[node] = parentvals[1] .* parentvals[2]
        else
            vals[node] = parentvals[1] * parentvals[2]
        end
    elseif (op == :lincomb)
        vals[node] = sum(graph.coeffs[node].*parentvals);
    elseif (op == :ldiv)
        if (parentsvals[1] isa Vector)
            vals[node] = parentvals[1] .\ parentvals[2]
        else
            vals[node] = parentvals[1] \ parentvals[2]
        end
    else
        error("Unknown operation")
    end
    return nothing
end

# Initiate the Dict vals with different types via dispatch. Create if needed
function init_vals_eval_graph!(graph, x, vals, input)
    if vals == nothing
        T_comp = promote_type(eltype(graph), typeof(x))
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I] = one(x)
    vals[input] = x
    return vals
end
function init_vals_eval_graph!(graph, x::AbstractVector, vals, input)
    if vals == nothing
        T_elem = promote_type(eltype(graph), eltype(x))
        if x isa Vector
            T_comp = Vector{T_elem}
        else
            T_comp = AbstractVector{T_elem}
        end
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I] = one.(x)
    vals[input] = x
    return vals
end
function init_vals_eval_graph!(graph, x::AbstractMatrix, vals, input)
    if vals == nothing
        T_elem = promote_type(eltype(graph), eltype(x))
        if x isa Matrix
            T_comp = Matrix{T_elem}
        elseif x isa SparseMatrixCSC
            T_comp = SparseMatrixCSC{T_elem,Int}
        else
            T_comp = AbstractMatrix{T_elem}
        end
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I] = one(x)
    vals[input] = x
    return vals
end

## Compute derivatives w.r.t. coefficients in the graph

# Computes Jacobian J = dZ(x_i)/dc, with respect to the coefficients
# c given in the vector cref, and points x_i in the vector x.
# If the objective funcion is the sum of squared residuals, then the gradient
# is sum_i (Z(x_i)-f(x_i)) dZ(x_i)/dc, i.e., J^T times vector of residual
# values r_i=Z(x_i)-f(x_i).
"""
    J=eval_jac(
        graph,
        x,
        cref;
        vals = nothing,
        input = :A,
        output = size(graph.outputs, 1),
    )

Computes Jacobian `J = dZ(x_i)/dc`, with respect to the coefficients given in
the vector `cref`, and points in `x`. See [`eval_graph`](@ref) for description
of `vals` and `output`.
"""
function eval_jac(
    graph,
    x,
    cref;
    vals = nothing,
    input = :A,
    output = size(graph.outputs, 1),
)
    if (cref isa Tuple) #Workaround for single coefficient
        cref = [cref]
    end

    # T = promote_type(eltype(x),eltype(graph))
    comporder = get_topo_order(graph, input = input)[1]
    if vals == nothing
        vals = init_vals_eval_graph!(graph, x, vals, input)
        eval_graph(
            graph,
            x,
            vals = vals,
            input = input,
            output = output,
            comporder = comporder,
        )
    end
    T = eltype(valtype(vals))

    jac = zeros(T, length(x), length(cref))
    for (i, c) in enumerate(cref)
        jac[:, i] = eval_der(graph, x, c, vals, comporder, output)
    end
    return jac
end

function eval_der(graph, x, c, vals, comporder, output)
    # Computes derivatives J = dZ(x_i)/dc, with respect to the coefficient c
    # and points x_i in the vector x.
    T = eltype(valtype(vals))
    der = Dict{Symbol,Vector{T}}()
    # Allocate return vector
    der[graph.outputs[output]] = zeros(T, length(x))
    # Base case. Derivative is parent node value
    der[c[1]] = vals[graph.parents[c[1]][c[2]]]

    # Loop over topological sorting
    # Only compute if parent has a nonzero derivative
    for node in comporder
        if haskey(der, graph.parents[node][1]) ||
           haskey(der, graph.parents[node][2])
            der[node] = zeros(T, length(x))
            for i = 1:2 # Two parents. One or both has a nonzeros derivative
                curr_parent = graph.parents[node][i]
                other_parent = graph.parents[node][(i == 1) ? 2 : 1]
                if haskey(der, curr_parent)
                    der_prop!(
                        graph,
                        der,
                        vals,
                        node,
                        curr_parent,
                        other_parent,
                        i,
                    )
                end
            end
        end
    end
    return der[graph.outputs[output]]
end

# der_prop! -- Derivative forward propagation
function der_prop!(graph, der, vals, node, curr_parent, other_parent, i)
    op = graph.operations[node]
    if (op == :mult)
        der[node][:] = der[node][:] + vals[other_parent] .* der[curr_parent]
    elseif (op == :ldiv)
        if (i == 1)
            der[node][:] =
                der[node][:] -
                (vals[other_parent] ./ (vals[curr_parent] .^ 2)) .*
                der[curr_parent]
        elseif (i == 2)
            der[node][:] = der[node][:] + der[curr_parent] ./ vals[other_parent]
        end
    elseif (op == :lincomb)
        α = graph.coeffs[node]
        der[node][:] = der[node][:] + α[i] * der[curr_parent]
    else
        error("Unknown operation")
    end
    return nothing
end

## Compute running error
# eval_runerr()
# Modifies the Dict relerrs if passed as a kwarg
# Returns the output specificed in graph.outputs[output]
"""
    err=eval_runerr(
        graph,
        x;
        vals = nothing,
        relerrs = nothing,
        input = :A,
        output = size(graph.outputs, 1),
        mode = :bounds, # Can be :bounds, :rand, :estimate
        add_relerr = eps(),
    )

Provides the running error of the graph evaluated in `x`. See
[`eval_graph`](@ref) for meaning of `vals` and `output`. The kwarg `relerrs` is
an anologous variable for the running error estimates in each node. The kwarg
`mode` can be `:bounds`, `:rand`, `:estimate`, specifying if the code should
make estimates or bounds, or random error within the bound.
"""
function eval_runerr(
    graph,
    x;
    vals = nothing,
    relerrs = nothing,
    input = :A,
    output = size(graph.outputs, 1),
    mode = :bounds, # Can be :bounds, :rand, :estimate
    add_relerr = eps(),
)
    T = promote_type(eltype(x), eltype(graph))
    comporder = get_topo_order(graph; input = input)[1]
    if vals == nothing
        vals = init_vals_eval_graph!(graph, x, vals, input)
        eval_graph(
            graph,
            x,
            vals = vals,
            input = input,
            output = output,
            comporder = comporder,
        )
    end
    relerrs = init_relerrs_eval_runerr!(relerrs, T, x, add_relerr, input)

    for node in comporder
        parentval1 = vals[graph.parents[node][1]]
        parentval2 = vals[graph.parents[node][2]]
        nodeval = vals[node]
        error_prop!(
            graph,
            relerrs,
            nodeval,
            parentval1,
            parentval2,
            node,
            mode,
            add_relerr,
        )
    end
    return relerrs[graph.outputs[output]]
end

# error_prop! -- Error propagation with different types via dispatch
function error_prop!(
    graph,
    relerrs,
    nodeval,
    parentval1,
    parentval2,
    node,
    mode,
    add_relerr,
)
    op = graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]
    # Generic type: Treat as "scalar"
    if (op == :mult || op == :ldiv)
        relerrs[node] = relerrs[parent1] + relerrs[parent2]
    elseif (op == :lincomb)
        if (mode == :rand || mode == :estimate)
            if (nodeval == 0)
                nodeval = eps(BigFloat) * 100 # Hack to avoid division by zero
            end
            relerrs[node] =
                (
                    parentval1 * relerrs[parent1] +
                    parentval2 * relerrs[parent2]
                ) / nodeval
        else
            absval = abs(nodeval)
            if (absval == 0)
                absval = eps(BigFloat) * 100 # Hack to avoid division by zero
            end
            relerrs[node] =
                (
                    abs(parentval1) * relerrs[parent1] +
                    abs(parentval2) * relerrs[parent2]
                ) / absval
        end
    else
        error("Unknown operation")
    end
    if (mode == :rand)
        relerrs[node] = (1 - 2 * rand()) * relerrs[node]
    end
    # Introduce round-off errors.
    relerrs[node] = relerrs[node] + add_relerr
    return nothing
end
function error_prop!(
    graph,
    relerrs,
    nodeval::AbstractVector,
    parentval1::AbstractVector,
    parentval2::AbstractVector,
    node,
    mode,
    add_relerr,
)
    # Vectors: Operate on all points in parallell
    op = graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]
    if (op == :mult || op == :ldiv)
        relerrs[node] = relerrs[parent1] + relerrs[parent2]
    elseif (op == :lincomb)
        if (mode == :rand || mode == :estimate)
            nodeval[nodeval.==0] .= eps(BigFloat) * 100 # Avoid division by zero
            relerrs[node] =
                (
                    parentval1 .* relerrs[parent1] +
                    parentval2 .* relerrs[parent2]
                ) ./ nodeval
        else
            absval = abs.(nodeval)
            absval[absval.==0] .= eps(BigFloat) * 100 # Avoid division by zero
            relerrs[node] =
                (
                    abs.(parentval1) .* relerrs[parent1] +
                    abs.(parentval2) .* relerrs[parent2]
                ) ./ absval
        end
    else
        error("Unknown operation")
    end
    if (mode == :rand)
        relerrs[node] = (1 .- 2 * rand(length(relerrs[node]))) .* relerrs[node]
    end
    # Introduce round-off errors.
    relerrs[node] = relerrs[node] .+ add_relerr
    return nothing
end

# Initiate the Dict relerr with different types via dispatch. Create if needed
function init_relerrs_eval_runerr!(
    relerrs,
    T,
    x,
    add_relerr,
    input
)
    if relerrs == nothing
        T_big = big(T)
        relerrs = Dict{Symbol,T_big}()
    end
    relerrs[:I] = add_relerr
    relerrs[input] = add_relerr
    return relerrs
end
function init_relerrs_eval_runerr!(
    relerrs,
    T,
    x::AbstractVector,
    add_relerr,
    input,
)
    if relerrs == nothing
        T_big = big(T)
        relerrs = Dict{Symbol,Vector{T_big}}()
    end
    relerrs[:I] = zeros(T_big, length(x))
    relerrs[:I][:] .= add_relerr
    relerrs[input] = zeros(T_big, length(x))
    relerrs[input][:] .= add_relerr
    return relerrs
end
