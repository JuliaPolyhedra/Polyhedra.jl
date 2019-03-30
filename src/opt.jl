export AbstractPolyhedraOptimizer

struct PolyhedraOptSet{T, RepT <: Rep{T}} <: MOI.AbstractVectorSet
    rep::RepT
end

function JuMP.build_constraint(error_fun::Function, func, set::Rep)
    return JuMP.build_constraint(error_fun, func, PolyhedraOptSet(set))
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
    optimizer.objective_func = sparsevec(indices, coefs)
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

"""
    isempty(p::Rep, solver::JuMP.OptimizerFactory=Polyhedra.solver(p))

Check whether the polyhedron `p` is empty by using the solver `solver`.
"""
function Base.isempty(p::Rep{T}, solver::Solver=Polyhedra.solver(p)) where {T}
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
    model = JuMP.Model(solver)
    x = JuMP.@variable(model, [1:fulldim(p)])
    JuMP.@constraint(model, x in p)
    JuMP.optimize!(model)
    term = termination_status(model)
    if term == MOI.OPTIMAL
        return false
    elseif term == MOI.INFEASIBLE
        return true
    else
        error("Cannot determine whether the polyhedron is empty or not because",
              " the linear program terminated with status $term.")
    end
end
