export AbstractPolyhedraOptimizer

struct PolyhedraOptSet{T, RepT <: Rep{T}} <: MOI.AbstractVectorSet
    rep::RepT
end
MOI.dimension(set::PolyhedraOptSet) = fulldim(set.rep)

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
function MOI.get(optimizer::AbstractPolyhedraOptimizer, ::MOI.ObjectiveSense)
    return optimizer.objective_sense
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
function MOI.delete(optimizer::AbstractPolyhedraOptimizer{T},
                    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, <:Union{MOI.EqualTo{T}, MOI.LessThan{T}}}) where T
    optimizer.feasible_set = nothing # Invalidated by the deleted constraint
    return MOI.delete(optimizer.lphrep, ci)
end
function MOI.modify(optimizer::AbstractPolyhedraOptimizer{T},
                    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, <:Union{MOI.EqualTo{T}, MOI.LessThan{T}}},
                    change::MOI.AbstractFunctionModification) where T
    optimizer.feasible_set = nothing # Invalidated by the deleted constraint
    return MOI.modify(optimizer.lphrep.model, ci, change)
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

const NO_SOLVER_HELP = """
To provide a solver to a polyhedron, first select a solver from https://www.juliaopt.org/JuMP.jl/stable/installation/#Getting-Solvers-1.
If you choose for instance `GLPK`, do `using GLPK; solver = GLPK.Optimizer`.
Then provide the solver to the library. For instance, with the default library, do `lib = DefaultLibrary{Float64}(solver)`
or if you use an external library, say `QHull`, do `lib = QHull.Library(solver)`.
Then when you create the polyhedron, say from a representation `rep`, do `polyhedron(rep, lib)`.
"""

# Inspired by `JuMP.set_optimizer`
function layered_optimizer(solver)
    solver === nothing && error("No solver specified.\n", NO_SOLVER_HELP)
    optimizer = MOI.instantiate(solver)
    T = coefficient_type(optimizer)
    if !MOIU.supports_default_copy_to(optimizer, false)
        universal_fallback = MOIU.UniversalFallback(_MOIModel{T}())
        optimizer = MOIU.CachingOptimizer(universal_fallback, optimizer)
    end
    optimizer = MOI.Bridges.full_bridge_optimizer(optimizer, T)
    MOI.Bridges.add_bridge(optimizer, PolyhedraToLPBridge{T})
    return optimizer, T
end

_optimize!(model::JuMP.Model) = JuMP.optimize!(model)
_optimize!(model::MOI.ModelLike) = MOI.optimize!(model)

function _unknown_status(model, status, message)
    error("Solver returned ", status, " when ", message, " Solver specific status: ", MOI.get(model, MOI.RawStatusString()))
end
function is_feasible(model, message)
    @assert MOI.get(model, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    _optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    # The objective sense is `MOI.FEASIBILITY_SENSE` so
    # `INFEASIBLE_OR_UNBOUNDED` means infeasible because it cannot be
    # unbounded
    if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        return false
    elseif status == MOI.OPTIMAL
        return true
    else
        _unknown_status(model, status, message)
    end
end

"""
    isempty(p::Rep, solver=Polyhedra.linear_objective_solver(p))

Check whether the polyhedron `p` is empty by using the solver `solver`.
"""
function Base.isempty(p::Rep, solver=Polyhedra.linear_objective_solver(p))
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
    x, cx = MOI.add_constrained_variables(model, PolyhedraOptSet(p))
    return !is_feasible(model, "trying to determine whether the polyhedron is empty.")
end
