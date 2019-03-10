using JuMP

"""
    hrep(model::JuMP.Model)

Builds an H-representation from the feasibility set of the JuMP model `model`.
Note that if non-linear constraint are present in the model, they are ignored.
"""
hrep(model::JuMP.Model) = LPHRepresentation(model)

# Returns affine constraint lower and upper bounds, all as dense vectors
function prepConstrBounds(m::Model)

    # Create dense affine constraint bound vectors
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # -Inf means no lower bound, +Inf means no upper bound
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    @inbounds for ind in 1:numRows
        rowlb[ind] = linconstr[ind].lb
        rowub[ind] = linconstr[ind].ub
    end

    return rowlb, rowub
end

# Converts all the affine constraints into a sparse column-wise
# matrix of coefficients.
function prepConstrMatrix(m::Model)

    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # Calculate the maximum number of nonzeros
    # The actual number may be less because of cancelling or
    # zero-coefficient terms
    nnz = 0
    for c in 1:numRows
        nnz += length(linconstr[c].terms.coeffs)
    end
    # Non-zero row indices
    I = Array{Int}(nnz)
    # Non-zero column indices
    J = Array{Int}(nnz)
    # Non-zero values
    V = Array{Float64}(nnz)

    # Fill it up!
    # Number of nonzeros seen so far
    nnz = 0
    for c in 1:numRows
        # Check that no coefficients are NaN/Inf
        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars   = linconstr[c].terms.vars
        # Check that variables belong to this model
        if !verify_ownership(m, vars)
            throw(VariableNotOwnedError("constraint"))
        end
        # Record all (i,j,v) triplets
        @inbounds for ind in 1:length(coeffs)
            nnz += 1
            I[nnz] = c
            J[nnz] = vars[ind].col
            V[nnz] = coeffs[ind]
        end
    end

    # sparse() handles merging duplicate terms and removing zeros
    A = sparse(I,J,V,numRows,m.numCols)
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
