using Test, LinearAlgebra, SparseArrays
using JuMP
using Polyhedra

mutable struct MockOptimizer{T} <: Polyhedra.AbstractPolyhedraOptimizer{T}
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

    optimize!::Function
    status::MOI.TerminationStatusCode
    solution::Union{AbstractVector{T}, Nothing}

    function MockOptimizer{T}(optimize!::Function) where T
        return new(LPHRep(Polyhedra._MOIModel{T}()), nothing, nothing,
                   MOI.FEASIBILITY_SENSE, nothing, zero(T),
                   optimize!, MOI.OPTIMIZE_NOT_CALLED, nothing)
    end
end
coefficient_type(::MockOptimizer{T}) where {T} = T

function MOI.empty!(mock::MockOptimizer{T}) where T
    mock.lphrep = LPHRep(Polyhedra._MOIModel{T}())
    mock.rep = nothing
    mock.feasible_set = nothing
    mock.objective_sense = MOI.FEASIBILITY_SENSE
    mock.objective_func = nothing
    mock.objective_constant = zero(T)
    mock.status = MOI.OPTIMIZE_NOT_CALLED
    mock.solution = nothing
end
function MOI.is_empty(mock::MockOptimizer{T}) where T
    MOI.is_empty(mock.lphrep.model) &&
    mock.rep === nothing &&
    mock.feasible_set === nothing &&
    mock.objective_sense == MOI.FEASIBILITY_SENSE &&
    mock.objective_func === nothing &&
    iszero(mock.objective_constant) &&
    mock.status == MOI.OPTIMIZE_NOT_CALLED &&
    mock.solution === nothing
end

function MOI.optimize!(mock::MockOptimizer)
    mock.optimize!(mock)
end

MOI.get(mock::MockOptimizer, ::MOI.TerminationStatus) = mock.status
function MOI.get(mock::MockOptimizer, ::MOI.ResultCount)
    if mock.status == MOI.OPTIMAL || mock.status == MOI.DUAL_INFEASIBLE
        return 1
    else
        return 0
    end
end
function MOI.get(mock::MockOptimizer{T}, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                         <:Union{MOI.EqualTo{T},
                                                 MOI.LessThan{T}}}) where T
    return MOIU.get_fallback(mock, attr, ci)
end
function MOI.get(mock::MockOptimizer, ::MOI.ObjectiveValue)
    if mock.status == MOI.OPTIMAL
        return mock.objective_func ⋅ mock.solution + mock.objective_constant
    elseif mock.status == MOI.DUAL_INFEASIBLE
        return mock.objective_func ⋅ mock.solution
    else
        error("No objective value available when termination status is $(mock.status).")
    end
end
function MOI.get(mock::MockOptimizer, ::MOI.PrimalStatus)
    if mock.status == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    elseif mock.status == MOI.DUAL_INFEASIBLE
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.NO_SOLUTION
    end
end
function MOI.get(mock::MockOptimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if mock.status != MOI.OPTIMAL && mock.status != MOI.DUAL_INFEASIBLE
        error("No primal value available when termination status is $(mock.status).")
    end
    return mock.solution[vi.value]
end

function polyhedra_model_test(factory::JuMP.OptimizerFactory)
    model = Model(factory)
    n = 3
    @variable(model, 0 <= x[1:n] <= 1)
    @variable(model, y)
    @constraint(model, 2y == 2)
    csum = @constraint(model, sum(x[i] for i=1:n) == 1)
    c12 = @constraint(model, x[1] - x[2] <= 1)
    c23 = @constraint(model, x[2] - x[3] >= 1)
    @objective(model, Max, sum(x) - 2y)

    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) == -1.0
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value.(x) ≈ [0.0, 1.0, 0.0]
    @test value(y) ≈ 1.0
    @test value(csum) ≈ 1.0
    @test value(c12) ≈ -1.0
    @test value(c23) ≈ 1.0

    # TODO dual
end

@testset "AbstractPolyhedraOptimizer" begin
    polyhedra_model_test(with_optimizer(
        MockOptimizer{Float64},
        mock -> begin
            mock.status = MOI.OPTIMAL
            mock.solution = [0.0, 1.0, 0.0, 1.0]
        end))
end
