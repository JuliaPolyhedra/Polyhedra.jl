using JuMP

function LPHRepresentation(model::JuMP.Model)
    # Inspired from Joey Huchette's code in ConvexHull.jl
    A = JuMP.prepConstrMatrix(model)
    c = JuMP.prepAffObjective(model)
    lb, ub = JuMP.prepConstrBounds(model)
    l, u = model.colLower, model.colUpper

    m, n = size(A)
    @assert m == length(lb) == length(ub)
    @assert model.nlpdata == nothing
    @assert isempty(model.quadconstr)
    @assert isempty(model.sosconstr)

    LPHRepresentation(A, l, u, lb, ub)
end

function polyhedron(model::JuMP.Model, lib::PolyhedraLibrary=getlibraryfor(model.numCols, Float64))
    polyhedron(LPHRepresentation(model), lib)
end
