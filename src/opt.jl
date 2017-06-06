using MathProgBase

export AbstractPolyhedraModel

abstract type AbstractPolyhedraModel <: MathProgBase.AbstractLinearQuadraticModel end

# see the cheat in lpqp_to_polyhedra
#function PolyhedraModel(solver::MathProgBase.AbstractMathProgSolver)
#  error("PolyhedraModel not implemented for solver $solver")
#end

#creates ambiguity with (::SimpleVRepMode, ::Vrep) when a polyhedron is given
#loadproblem!(m::AbstractPolyhedraModel, hrep::HRep, c, sense) = error("loadproblem! not implemented")

struct LinprogSolution
    status
    objval
    sol
    attrs
end

function MathProgBase.linprog{N}(c::AbstractVector, p::Rep{N}, solver::MathProgBase.AbstractMathProgSolver = defaultLPsolverfor(p))
    m = PolyhedraModel(solver)
    if N != length(c)
        println("length of objective does not match dimension of polyhedron")
    end
    loadproblem!(m, p, c, :Min)
    optimize!(m)
    stat = status(m)
    if stat == :Optimal
        return LinprogSolution(stat, getobjval(m), getsolution(m), Dict())
    elseif stat == :Unbounded
        attrs = Dict()
        try
            attrs[:unboundedray] = getunboundedray(m)
        catch
            warn("Problem is unbounded, but unbounded ray is unavailable; check that the proper solver options are set.")
        end
        return LinprogSolution(stat, nothing, [], attrs)
    else
        return LinprogSolution(stat, nothing, [], Dict())
    end
end

function Base.isempty{N,T}(p::Polyhedron{N,T}, solver::MathProgBase.AbstractMathProgSolver = defaultLPsolverfor(p))
    linprog(zeros(T, N), p, solver).status == :Infeasible
end
