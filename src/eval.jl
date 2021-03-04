
export eval_graph,eval_jac,eval_runerr,adjust_for_errtype!

# Evaluation of graph, Jacobian, and running error
# eval_graph()
# eval_jac()
# adjust_for_errtype!()
# eval_runerr()

# init_vals_eval_graph!()
# init_relerrs_eval_runerr!()

## Evaluating the graph

# eval_graph()
# Modifies the Dict vals if passed as a kwarg
# Returns the output specificed in graph.outputs[output]
function eval_graph(graph, x; vals=nothing, output=size(graph.outputs,1), comporder=nothing)
    vals = init_vals_eval_graph!(graph, x, vals)
    if comporder == nothing
        comporder = get_computation_order(graph)[1]
    end
    for node = comporder
        parentval1=vals[graph.parents[node][1]]
        parentval2=vals[graph.parents[node][2]]
        carry_out!(graph,vals,parentval1,parentval2,node)
    end
    return vals[graph.outputs[output]]
end

# carry_out! -- Perform the operations with different types via dispatch
function carry_out!(graph, vals, parentval1, parentval2, node)
    # Generic type: Treat as "scalar"
    op=graph.operations[node]
    if (op==:mult)
        vals[node]=parentval1*parentval2
    elseif (op==:lincomb)
        (α1,α2)=graph.coeffs[node]
        vals[node]=α1*parentval1+α2*parentval2
    elseif (op==:ldiv)
        vals[node]=parentval1\parentval2
    else
        error("Unknown operation")
    end
    return nothing
end
function carry_out!(graph, vals, parentval1::AbstractVector, parentval2::AbstractVector, node)
    # Vectors: Operate on all points in parallell
    op=graph.operations[node]
    if (op==:mult)
        vals[node]=parentval1 .* parentval2
    elseif (op==:lincomb)
        (α1,α2)=graph.coeffs[node]
        vals[node]=α1*parentval1+α2*parentval2
    elseif (op==:ldiv)
        vals[node]=parentval1 .\ parentval2
    else
        error("Unknown operation")
    end
    return nothing
end
function carry_out!(graph, vals, parentval1::AbstractMatrix, parentval2::AbstractMatrix, node)
    # Matrices: Optimized operations
    op=graph.operations[node]
    if (op==:mult)
        # result[:]=parentval1*parentval2
        vals[node]=similar(parentval2)
        mul!(vals[node],parentval1,parentval2)
        #BLAS.gemm!('N','N',1.0,parentval1,parentval2,0.0,result)
    elseif (op==:lincomb)
        (α1,α2)=graph.coeffs[node]
        vals[node]=copy(parentval2) # Assumes second argument is copy(parentval2)
        BLAS.axpby!(α1,parentval1,α2,vals[node])
        #result[:]=α1*parentval1+α2*parentval2 # equivalent
    elseif (op==:ldiv)
        vals[node]=parentval1\parentval2
    else
        error("Unknown operation")
    end
    return nothing
end

# Initiate the Dict vals with different types via dispatch. Create if needed
function init_vals_eval_graph!(graph, x, vals)
    if vals == nothing
        T_comp=promote_type(eltype(graph),typeof(x))
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I]=one(x)
    vals[:A]=x
    return vals
end
function init_vals_eval_graph!(graph, x::AbstractVector, vals)
    if vals == nothing
        T_elem=promote_type(eltype(graph),eltype(x))
        if x isa Vector
            T_comp=Vector{T_elem}
        else
            T_comp=AbstractVector{T_elem}
        end
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I]=one.(x)
    vals[:A]=x
    return vals
end
function init_vals_eval_graph!(graph, x::AbstractMatrix, vals)
    if vals == nothing
        T_elem=promote_type(eltype(graph),eltype(x))
        if x isa Matrix
            T_comp=Matrix{T_elem}
        elseif x isa SparseMatrixCSC
            T_comp=SparseMatrixCSC{T_elem,Int}
        else
            T_comp=AbstractMatrix{T_elem}
        end
        vals = Dict{Symbol,T_comp}()
    end
    vals[:I]=one(x)
    vals[:A]=x
    return vals
end


## Compute derifvatives w.r.t. coefficients in the graph

# Computes Jacobian J = dZ(x_i)/dc, with respect to the coefficients
# c given in the vector cref, and points x_i in the vector x.
# If the objective funcion is the sum of squared residuals, then the gradient
# is sum_i (Z(x_i)-f(x_i)) dZ(x_i)/dc, i.e., J^T times vector of residual
# values r_i=Z(x_i)-f(x_i).
function eval_jac(graph, x, cref; vals=nothing, output=size(graph.outputs,1))
    if (cref isa Tuple) #Workaround for single coefficient
        cref=[cref]
    end

    # T = promote_type(eltype(x),eltype(graph))
    comporder = get_computation_order(graph)[1]
    if vals==nothing
        vals = init_vals_eval_graph!(graph, x, vals)
        eval_graph(graph, x, vals=vals, output=output, comporder=comporder)
    end
    T = eltype(valtype(vals))


    jac = zeros(T,length(x),length(cref))
    for (i,c) = enumerate(cref)
        jac[:,i] = eval_der(graph, x, c, vals, comporder, output)
    end
    return jac
end

function eval_der(graph, x, c, vals, comporder, output)
    # Computes derivatives J = dZ(x_i)/dc, with respect to the coefficient c
    # and points x_i in the vector x.
    T = eltype(valtype(vals))
    der = Dict{Symbol,Vector{T}}()
    # Allocate return vector
    der[graph.outputs[output]] = zeros(T,length(x))
    # Base case. Derivative is parent node value
    der[c[1]] = vals[ graph.parents[c[1]][c[2]] ]

    # Loop over topological sorting
    # Only compute if parent has a nonzero derivative
    for node in comporder
        if haskey(der,graph.parents[node][1]) || haskey(der,graph.parents[node][2])
            der[node] = zeros(T,length(x))
            for i = 1:2 # Two parents. One or both has a nonzeros derivative
                curr_parent = graph.parents[node][i]
                other_parent = graph.parents[node][(i==1) ? 2 : 1]
                if haskey(der,curr_parent)
                    der_prop!(graph, der, vals, node, curr_parent, other_parent, i)
                end
            end
        end
    end
    return der[graph.outputs[output]]
end

# der_prop! -- Derivative forward propagation
function der_prop!(graph, der, vals, node, curr_parent, other_parent, i)
    op = graph.operations[node]
    if (op==:mult)
        der[node][:] = der[node][:] + vals[other_parent] .* der[curr_parent]
    elseif (op==:ldiv)
        if (i == 1)
            der[node][:] = der[node][:] - (vals[other_parent]./(vals[curr_parent].^2)) .* der[curr_parent]
       elseif (i==2)
           der[node][:] = der[node][:] + der[curr_parent]./vals[other_parent]
       end
    elseif (op==:lincomb)
        α=graph.coeffs[node]
        der[node][:] = der[node][:] + α[i]*der[curr_parent]
    else
        error("Unknown operation")
    end
    return nothing
end


# adjust_for_errtype!()
# JR can be Jacobian or vector of residuals.
# f is the objective function
# errtype can be :abserr or :relerr, although the former is doing nothing
function adjust_for_errtype!(JR, pts, f, errtype)
    # Assumes a Jacobian and residual-vector that is computed in absolute error.
    # Adjusts Jacobian and residual-vector from absolute, to relative error, i.e.,
    # objective function 1/2 sum_i (Z(x_i)-f(x_i))^2 / f(x_i)^2 and
    if (errtype==:abserr)
        # Do nothing
    elseif (errtype==:relerr)
        F = f.(pts)
        F[F.==0] .= eps()*100 # Hack to avoid division by zero
        JR[:] = Diagonal(1 ./ F)*JR
    else
        error("Unknown errtype '", errtype, "'.")
    end
    return JR
end

## Compute running error
# eval_runerr()
# Modifies the Dict relerrs if passed as a kwarg
# Returns the output specificed in graph.outputs[output]
function eval_runerr(graph, x; vals=nothing, relerrs=nothing,
                    output=size(graph.outputs,1),
                    mode=:bounds, # Can be :bounds, :rand, :estimate
                    add_relerr=eps())
    T = promote_type(eltype(x),eltype(graph))
    comporder = get_computation_order(graph)[1]
    if vals==nothing
        vals = init_vals_eval_graph!(graph, x, vals)
        eval_graph(graph, x, vals=vals, output=output, comporder=comporder)
    end
    relerrs = init_relerrs_eval_runerr!(relerrs, T, x, add_relerr)

    for node = comporder
        parentval1=vals[graph.parents[node][1]]
        parentval2=vals[graph.parents[node][2]]
        nodeval=vals[node]
        error_prop!(graph, relerrs, nodeval, parentval1, parentval2, node, mode, add_relerr)
    end
    return relerrs[graph.outputs[output]]
end

# error_prop! -- Error propagation with different types via dispatch
function error_prop!(graph, relerrs, nodeval, parentval1, parentval2, node, mode, add_relerr)
    op=graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]
    # Generic type: Treat as "scalar"
    if (op==:mult || op==:ldiv)
        relerrs[node]=relerrs[parent1]+relerrs[parent2]
    elseif (op==:lincomb)
        if (mode == :rand || mode == :estimate)
            if (nodeval==0)
                nodeval=eps(BigFloat)*100 # Hack to avoid division by zero
            end
            relerrs[node]=(parentval1*relerrs[parent1]
                        + parentval2*relerrs[parent2])/nodeval
        else
            absval=abs(nodeval)
            if (absval==0)
                absval=eps(BigFloat)*100 # Hack to avoid division by zero
            end
            relerrs[node]=(abs(parentval1)*relerrs[parent1]
                        + abs(parentval2)*relerrs[parent2])/absval
        end
    else
        error("Unknown operation")
    end
    if (mode == :rand)
        relerrs[node]=(1-2*rand())*relerrs[node]
    end
    # Introduce round-off errors.
    relerrs[node] = relerrs[node] + add_relerr
    return nothing
end
function error_prop!(graph, relerrs, nodeval::AbstractVector, parentval1::AbstractVector, parentval2::AbstractVector, node, mode, add_relerr)
    # Vectors: Operate on all points in parallell
    op=graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]
    if (op==:mult || op==:ldiv)
        relerrs[node]=relerrs[parent1]+relerrs[parent2]
    elseif (op==:lincomb)
        if (mode == :rand || mode == :estimate)
            nodeval[ nodeval.==0 ] .= eps(BigFloat)*100 # Hack to avoid division by zero
            relerrs[node]=(parentval1 .* relerrs[parent1]
                        + parentval2 .* relerrs[parent2]) ./ nodeval
        else
            absval=abs.(nodeval)
            absval[ absval.==0 ] .= eps(BigFloat)*100 # Hack to avoid division by zero
            relerrs[node]=(abs.(parentval1) .* relerrs[parent1]
                        + abs.(parentval2) .* relerrs[parent2]) ./ absval
        end
    else
        error("Unknown operation")
    end
    if (mode == :rand)
        relerrs[node]=(1 .- 2*rand(length(relerrs[node]))) .* relerrs[node]
    end
    # Introduce round-off errors.
    relerrs[node] = relerrs[node] .+ add_relerr
    return nothing
end
function error_prop!(graph, relerrs, nodeval::AbstractMatrix, parentval1::AbstractMatrix, parentval2::AbstractMatrix, node, mode, add_relerr)
    # Matrices
    op=graph.operations[node]
    parent1 = graph.parents[node][1]
    parent2 = graph.parents[node][2]
    if (op==:mult) || (op==:ldiv)
        relerrs[node]=NaN # TODO: The formula is more complicated
    elseif (op==:lincomb)
        if (mode == :rand || mode == :estimate)
            relerrs[node]=NaN # TODO: A signed scalar value representing the error
        else
            valnorm=norm(nodeval)
            if (valnorm==0)
                valnorm=eps(BigFloat)*100 # Hack to avoid division by zero
            end
            relerrs[node]=(norm(parentval1)*relerrs[parent1]
                        + norm(parentval2)*relerrs[parent2])/valnorm
        end
    else
        error("Unknown operation")
    end
    if (mode == :rand)
        relerrs[node]=(1-2*rand())*relerrs[node] #NOTE: Okey like this if scalars
    end
    # Introduce round-off errors.
    relerrs[node] = relerrs[node] + add_relerr
    return nothing
end

# Initiate the Dict relerr with different types via dispatch. Create if needed
function init_relerrs_eval_runerr!(relerrs, T, x, add_relerr)
    if relerrs==nothing
        T_big=big(T)
        relerrs=Dict{Symbol,T_big}()
    end
    relerrs[:I]=add_relerr
    relerrs[:A]=add_relerr
    return relerrs
end
function init_relerrs_eval_runerr!(relerrs, T, x::AbstractVector, add_relerr)
    if relerrs==nothing
        T_big=big(T)
        relerrs=Dict{Symbol,Vector{T_big}}()
    end
    relerrs[:I]=zeros(T_big,length(x))
    relerrs[:I][:] .= add_relerr
    relerrs[:A]=zeros(T_big,length(x))
    relerrs[:A][:] .= add_relerr
    return relerrs
end



;
#
