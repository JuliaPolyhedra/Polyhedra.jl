export AbstractPolyhedraOptimizer

abstract type AbstractPolyhedraOptimizer{T} <: MOI.AbstractOptimizer end

function MOI.copy_to(dest::AbstractPolyhedraOptimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(dest, src; kws...)
end
MOI.Utilities.supports_default_copy_to(optimizer::AbstractPolyhedraOptimizer, copy_names::Bool) = true

struct PolyhedraOptSet{T, RepT <: Rep{T}} <: MOI.AbstractVectorSet
    rep::RepT
end

function MOI.add_variable(optimizer::AbstractPolyhedraOptimizer)
    return MOI.add_variable(optimizer.lphrep.model)
end
function MOI.add_variables(optimizer::AbstractPolyhedraOptimizer, n)
    return MOI.add_variables(optimizer.lphrep.model, n)
end

function MOI.add_constraint(optimizer::AbstractPolyhedraOptimizer{T},
                            func::MOI.ScalarAffineFunction{T},
                            set::Union{MOI.EqualTo{T}, MOI.LessThan{T}})
    return MOI.add_constraint(optimizer.lphrep.model, func, set)
end
function MOI.add_constraint(optimizer::AbstractPolyhedraOptimizer,
                            vov::MOI.VectorOfVariables, rep::PolyhedraOptSet)
    if vov.variables != MOI.get(optimizer.lphrep, MOI.ListOfVariableIndices())
        error("Cannot only add VectorOfVariables polyhedra constraints with all variables in creation order.")
    end
    if optimizer.rep !== nothing
        error("Cannot only add one polyhedra constraint.")
    end
    optimizer.rep = rep
    return MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(rep)}(1)
end

# see the cheat in lpqp_to_polyhedra
#function PolyhedraModel(solver::Solver)
#  error("PolyhedraModel not implemented for solver $solver")
#end

#creates ambiguity with (::VRepMode, ::Vrep) when a polyhedron is given
#loadproblem!(m::AbstractPolyhedraModel, hrep::HRep, c, sense) = error("loadproblem! not implemented")

struct LinprogSolution
    status
    objval
    sol
    attrs
end

"""
    linprog(c::AbstractVector, p::Rep, solver::MathProgBase.AbstractMathProgSolver=Polyhedra.solver(p))

Solve the minimization of the objective ``\\langle c, x \\rangle`` over the polyhedron `p`.
"""
function MPB.linprog(c::AbstractVector, p::Rep, solver::Solver=Polyhedra.solver(p))
    m = PolyhedraModel(solver)
    if fulldim(p) != length(c)
        throw(DimensionMismatch("length of objective does not match dimension of polyhedron"))
    end
    MPBSI.loadproblem!(m, p, c, :Min)
    MPBSI.optimize!(m)
    stat = MPBSI.status(m)
    if stat == :Optimal
        return LinprogSolution(stat, MPBSI.getobjval(m), MPBSI.getsolution(m), Dict())
    elseif stat == :Unbounded
        attrs = Dict()
        try
            attrs[:unboundedray] = MPBSI.getunboundedray(m)
        catch
            @warn "Problem is unbounded, but unbounded ray is unavailable; check that the proper solver options are set."
        end
        return LinprogSolution(stat, nothing, [], attrs)
    else
        return LinprogSolution(stat, nothing, [], Dict())
    end
end

"""
    isempty(p::Rep, solver::MathProgBase.AbstractMathProgSolver=Polyhedra.solver(p))

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
    return MPB.linprog(zeros(T, fulldim(p)), p, solver).status == :Infeasible
end
