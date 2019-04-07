export AbstractPolyhedraOptimizer

struct PolyhedraOptSet{T, RepT <: Rep{T}} <: MOI.AbstractVectorSet
    rep::RepT
end

function JuMP.build_constraint(error_fun::Function, func, set::Rep)
    return JuMP.BridgeableConstraint(
        JuMP.build_constraint(error_fun, func, PolyhedraOptSet(set)),
        PolyhedraToLPBridge)
end

abstract type AbstractPolyhedraOptimizer{T} <: MOI.AbstractOptimizer end

function MOI.copy_to(dest::AbstractPolyhedraOptimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(dest, src; kws...)
end
MOI.Utilities.supports_default_copy_to(optimizer::AbstractPolyhedraOptimizer, copy_names::Bool) = true

function MOI.add_variable(optimizer::AbstractPolyhedraOptimizer)
    return MOI.add_variable(optimizer.lphrep.model)
end
function MOI.add_variables(optimizer::AbstractPolyhedraOptimizer, n)
    return MOI.add_variables(optimizer.lphrep.model, n)
end

function MOI.supports(::AbstractPolyhedraOptimizer{T},
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.SingleVariable},
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}}) where T
    return true
end
function MOI.set(optimizer::AbstractPolyhedraOptimizer, ::MOI.ObjectiveSense,
                 sense::MOI.OptimizationSense)
    optimizer.objective_sense = sense
end
function MOI.set(optimizer::AbstractPolyhedraOptimizer{T}, ::MOI.ObjectiveFunction,
                 func::MOI.SingleVariable) where T
    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
            convert(MOI.ScalarAffineFunction{T}, func))
end
function MOI.set(optimizer::AbstractPolyhedraOptimizer, ::MOI.ObjectiveFunction,
                 func::MOI.ScalarAffineFunction)
    indices = [term.variable_index.value for term in func.terms]
    coefs = [term.coefficient for term in func.terms]
    optimizer.objective_func = sparsevec(indices, coefs, fulldim(optimizer.lphrep))
    optimizer.objective_constant = func.constant
end

function MOI.supports_constraint(::AbstractPolyhedraOptimizer{T},
                                 ::Type{MOI.ScalarAffineFunction{T}},
                                 ::Type{<:Union{MOI.EqualTo{T}, MOI.LessThan{T}}}) where T
    return true
end
function MOI.get(optimizer::AbstractPolyhedraOptimizer{T},
                 attr::Union{MOI.ConstraintFunction,
                             MOI.ConstraintSet},
                 ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                         <:Union{MOI.EqualTo{T},
                                                 MOI.LessThan{T}}}) where T
    return MOI.get(optimizer.lphrep.model, attr, ci)
end
function MOI.add_constraint(optimizer::AbstractPolyhedraOptimizer{T},
                            func::MOI.ScalarAffineFunction{T},
                            set::Union{MOI.EqualTo{T}, MOI.LessThan{T}}) where T
    optimizer.feasible_set = nothing # Invalidated by the new constraint
    return MOI.add_constraint(optimizer.lphrep, func, set)
end
function MOI.supports_constraint(::AbstractPolyhedraOptimizer,
                                 ::Type{MOI.VectorOfVariables},
                                 ::Type{<:PolyhedraOptSet})
    return true
end
function MOI.add_constraint(optimizer::AbstractPolyhedraOptimizer,
                            vov::MOI.VectorOfVariables, rep::PolyhedraOptSet)
    if vov.variables != MOI.get(optimizer.lphrep.model, MOI.ListOfVariableIndices())
        error("Cannot only add VectorOfVariables polyhedra constraints with all variables in creation order.")
    end
    if optimizer.rep !== nothing
        error("Cannot only add one polyhedra constraint.")
    end
    optimizer.feasible_set = nothing # Invalidated by the new constraint
    optimizer.rep = rep.rep
    return MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(rep)}(1)
end

coefficient_type(::MOI.ModelLike) = Float64

# Inspired by `JuMP.set_optimizer`
function layered_optimizer(factory::JuMP.OptimizerFactory)
    optimizer = factory()
    T = coefficient_type(optimizer)
    if !MOIU.supports_default_copy_to(optimizer, false)
        universal_fallback = MOIU.UniversalFallback(_MOIModel{T}())
        optimizer = MOIU.CachingOptimizer(universal_fallback, optimizer)
    end
    optimizer = MOI.Bridges.full_bridge_optimizer(optimizer, T)
    MOI.Bridges.add_bridge(optimizer, PolyhedraToLPBridge{T})
    return optimizer, T
end

"""
    isempty(p::Rep, solver::JuMP.OptimizerFactory=Polyhedra.linear_objective_solver(p))

Check whether the polyhedron `p` is empty by using the solver `solver`.
"""
function Base.isempty(p::Rep, solver::Solver=Polyhedra.linear_objective_solver(p))
    N = fulldim(p)
    if N == -1
        if p isa VRepresentation || (p isa Polyhedron && fulldim(vrep(p)) == -1)
            # Empty V-representation means empty polyhedron
            return true
        else
            # Empty H-representation means the polyhedron the full space
            return false
        end
    end
    model, T = layered_optimizer(solver)
    x = MOI.add_variables(model, fulldim(p))
    MOI.add_constraint(model, MOI.VectorOfVariables(x), PolyhedraOptSet(p))
    MOI.optimize!(model)
    term = MOI.get(model, MOI.TerminationStatus())
    if term == MOI.OPTIMAL
        return false
    elseif term == MOI.INFEASIBLE || term == MOI.INFEASIBLE_OR_UNBOUNDED
        # The problem has no objective so it cannot be unbounded so
        # the `MOI.INFEASIBLE_OR_UNBOUNDED` status is accepted
        return true
    else
        error("Cannot determine whether the polyhedron is empty or not because",
              " the linear program terminated with status $term.")
    end
end
