using JuMP

"""
    hrep(model::JuMP.Model)

Builds an H-representation from the feasibility set of the JuMP model `model`.
Note that if non-linear constraint are present in the model, they are ignored.
"""
hrep(model::JuMP.Model) = LPHRepresentation(model)

function _get_constraint_ref_array(model::Model)
    ref_arr = vcat([
        all_constraints(model, c...)
        for c in list_of_constraint_types(model)
    ]...)
end

function _get_constraint_array(model::Model)
    constraint_arr = constraint_object.(
        _get_constraint_ref_array(model)
    )
end

# Returns affine constraint lower and upper bounds, all as dense vectors
function prepConstrBounds(m::Model)

    # Create dense affine constraint bound vectors
    linconstr = _get_constraint_array(m)
    numRows = length(linconstr)
    # -Inf means no lower bound, +Inf means no upper bound
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    @inbounds for ind in 1:numRows
        set = linconstr[ind].set
        if set isa MOI.GreaterThan
            rowlb[ind] = set.lower
        elseif set isa MOI.LessThan
            rowlb[ind] = set.upper
        end
    end

    return rowlb, rowub
end

# Converts all the affine constraints into a sparse column-wise
# matrix of coefficients.
function prepConstrMatrix(m::Model)

    linconstr = _get_constraint_array(m)
    # Number of constraints
    numRows = length(linconstr)
    # Number of variables
    numCols = num_variables(m)

    # Calculate the maximum number of nonzeros
    # The actual number may be less because of cancelling or
    # zero-coefficient terms
    nnz = 0
    for c in 1:numRows
        func = linconstr[c].func
        if func isa VariableRef
            nnz += 1
        elseif func isa GenericAffExpr
            nnz += length(func.terms)
        end
    end
    # Non-zero row indices
    I = Array{Int}(undef, nnz)
    # Non-zero column indices
    J = Array{Int}(undef, nnz)
    # Non-zero values
    V = Array{Float64}(undef, nnz)

    # Fill it up!
    # Number of nonzeros seen so far
    nnz = 0
    for c in 1:numRows
        func = linconstr[c].func
        if func isa VariableRef
            nnz += 1
            I[nnz] = c
            J[nnz] = func.index.value
            V[nnz] = 1.0
        elseif func isa GenericAffExpr
            vars, coeffs = zip(linconstr[c].func.terms...)
            @assert all(isfinite.(coeffs))
            # Record all (i,j,v) triplets
            @inbounds for ind in 1:length(coeffs)
                nnz += 1
                I[nnz] = c
                J[nnz] = vars[ind].index.value
                V[nnz] = coeffs[ind]
            end
        end
    end

    # sparse() handles merging duplicate terms and removing zeros
    A = sparse(I, J, V, numRows, numCols)
end

function LPHRepresentation(model::JuMP.Model)
    # Inspired from Joey Huchette's code in ConvexHull.jl
    A = prepConstrMatrix(model)
    lb, ub = prepConstrBounds(model)
    l, u = model.colLower, model.colUpper

    m, n = size(A)
    @assert m == length(lb) == length(ub)
    @assert model.nlpdata == nothing
    @assert isempty(model.quadconstr)
    @assert isempty(model.sosconstr)

    LPHRepresentation(A, l, u, lb, ub)
end

function polyhedron(model::JuMP.Model, lib::Library=default_library(FullDim{model.numCols}(), Float64))
    polyhedron(LPHRepresentation(model), lib)
end
