export VRepOptimizer

"""
    VRepOptimizer{T} <: AbstractPolyhedraOptimizer{T}

Linear Programming solver using the V-representation of the feasible set to find
the optimal solution.
"""
mutable struct VRepOptimizer{T} <: AbstractPolyhedraOptimizer{T}
    library::Union{Nothing, Library}
    lphrep::LPHRep{T}
    # Feasible set defined by user
    rep::Union{Rep{T}, Nothing}
    # Either `rep` or `polyhedron(lphrep, library)`.
    # It is kept between consecutive solve if not modified,
    # e.g. if only the objective is changed.
    feasible_set::Union{Rep{T}, Nothing}

    objective_sense::MOI.OptimizationSense
    objective_func::Union{SparseVector{T, Int64}, Nothing}
    objective_constant::T

    status::MOI.TerminationStatusCode
    solution::Union{AbstractVector{T}, Nothing}

    function VRepOptimizer{T}(library::Union{Nothing, Library} = nothing) where T
        new(library, LPHRep(_MOIModel{T}()), nothing, nothing,
            MOI.FEASIBILITY_SENSE, nothing, zero(T),
            MOI.OPTIMIZE_NOT_CALLED, nothing)
    end
end

coefficient_type(::VRepOptimizer{T}) where {T} = T
MOI.get(::VRepOptimizer, ::MOI.SolverName) = "VRep"

function MOI.empty!(lpm::VRepOptimizer{T}) where T
    lpm.lphrep = LPHRep(_MOIModel{T}())
    lpm.rep = nothing
    lpm.feasible_set = nothing
    lpm.objective_sense = MOI.FEASIBILITY_SENSE
    lpm.objective_func = nothing
    lpm.objective_constant = zero(T)
    lpm.status = MOI.OPTIMIZE_NOT_CALLED
    lpm.solution = nothing
end
function MOI.is_empty(lpm::VRepOptimizer{T}) where T
    MOI.is_empty(lpm.lphrep.model) &&
    lpm.rep === nothing &&
    lpm.feasible_set === nothing &&
    lpm.objective_sense == MOI.FEASIBILITY_SENSE &&
    lpm.objective_func === nothing &&
    iszero(lpm.objective_constant) &&
    lpm.status == MOI.OPTIMIZE_NOT_CALLED &&
    lpm.solution === nothing
end

function MOI.optimize!(lpm::VRepOptimizer{T}) where T
    if lpm.rep === nothing
        lpm.feasible_set = lpm.lphrep
    else
        if hasallhalfspaces(lpm.lphrep)
            error("Cannot provide both a polyhedral feasible set and additional constraints.")
        end
        lpm.feasible_set = lpm.rep
    end
    if lpm.feasible_set isa HRepresentation
        if lpm.library === nothing
            lpm.feasible_set = polyhedron(lpm.feasible_set)
        else
            lpm.feasible_set = polyhedron(lpm.feasible_set, lpm.library)
        end
    end
    if lpm.feasible_set isa VRepresentation
        prob = lpm.feasible_set
    else
        @assert lpm.feasible_set isa Polyhedron
        prob = vrep(lpm.feasible_set)
    end
    N = fulldim(prob)
    if !haspoints(prob) && !haslines(prob) && !hasrays(prob)
        lpm.status = MOI.INFEASIBLE
        lpm.solution = nothing
    elseif lpm.objective_sense == MOI.FEASIBILITY_SENSE
        lpm.status = MOI.OPTIMAL
        lpm.solution = first(points(prob))
    else
        better(a, b) = (lpm.objective_sense == MOI.MAX_SENSE ? a > b : a < b)
        _better(a, b) = (lpm.objective_sense == MOI.MAX_SENSE ? _gt(a, b) : _lt(a, b))
        bestobjval = zero(T)
        lpm.solution = nothing
        for r in allrays(prob)
            objval = lpm.objective_func ⋅ r
            if _better(objval, bestobjval)
                bestobjval = objval
                lpm.solution = coord(r)
            end
        end
        if lpm.solution !== nothing
            lpm.status = MOI.DUAL_INFEASIBLE
        else
            for p in points(prob)
                objval = lpm.objective_func ⋅ p
                if lpm.solution === nothing || better(objval, bestobjval)
                    bestobjval = objval
                    lpm.solution = p
                end
            end
            lpm.status = MOI.OPTIMAL
        end
        @assert lpm.solution !== nothing
    end
end

MOI.get(lpm::VRepOptimizer, ::MOI.TerminationStatus) = lpm.status
function MOI.get(lpm::VRepOptimizer, ::MOI.ResultCount)
    if lpm.status == MOI.OPTIMAL || lpm.status == MOI.DUAL_INFEASIBLE
        return 1
    else
        return 0
    end
end
function MOI.get(lpm::VRepOptimizer{T}, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                         <:Union{MOI.EqualTo{T},
                                                 MOI.LessThan{T}}}) where T
    return MOIU.get_fallback(lpm, attr, ci)
end
function MOI.get(lpm::VRepOptimizer, ::MOI.ObjectiveValue)
    if lpm.status == MOI.OPTIMAL
        return lpm.objective_func ⋅ lpm.solution + lpm.objective_constant
    elseif lpm.status == MOI.DUAL_INFEASIBLE
        return lpm.objective_func ⋅ lpm.solution
    else
        error("No objective value available when termination status is $(lpm.status).")
    end
end
function MOI.get(lpm::VRepOptimizer, ::MOI.PrimalStatus)
    if lpm.status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif lpm.status == MOI.DUAL_INFEASIBLE
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.NO_SOLUTION
    end
end
function MOI.get(lpm::VRepOptimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if lpm.status != MOI.OPTIMAL && lpm.status != MOI.DUAL_INFEASIBLE
        error("No primal value available when termination status is $(lpm.status).")
    end
    return lpm.solution[vi.value]
end
