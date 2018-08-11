export VRepSolver, VRepPolyhedraModel

"""
    VRepSolver

Linear Programming solver using the V-representation of the feasible set to find the optimal solution.
"""
struct VRepSolver <: MPB.AbstractMathProgSolver
end

mutable struct VRepPolyhedraModel <: AbstractPolyhedraModel
    vrep::Nullable{VRep}
    obj::Nullable{Vector}
    sense::Symbol

    objval::Nullable
    solution::Nullable{Vector}
    status::Symbol

    function VRepPolyhedraModel()
        new(nothing, nothing, :Feas, nothing, nothing, :Undecided)
    end
end

PolyhedraModel(::VRepSolver) = VRepPolyhedraModel()
MPBSI.LinearQuadraticModel(s::VRepSolver) = PolyhedraToLPQPBridge(PolyhedraModel(s))

function MPBSI.loadproblem!(lpm::VRepPolyhedraModel, vrep::VRep, obj, sense)
    if !(sense in [:Max, :Min])
        error("sense should be :Max or :Min")
    end
    lpm.vrep = vrep
    if isapproxzero(obj)
        lpm.sense = :Feas
    else
        lpm.obj = copy(obj)
        lpm.sense = sense
    end
end

function MPBSI.optimize!(lpm::VRepPolyhedraModel)
    if isnull(lpm.vrep)
        error("No problem loaded")
    end
    prob = get(lpm.vrep)
    N = fulldim(prob)
    T = MultivariatePolynomials.coefficienttype(prob)
    lpm.status = :Undecided
    if !haspoints(prob) && !haslines(prob) && !hasrays(prob)
        lpm.status = :Infeasible
    else
        if lpm.sense in [:Max, :Min]
            better(a, b) = (lpm.sense == :Max ? a > b : a < b)
            _better(a, b) = (lpm.sense == :Max ? _gt(a, b) : _lt(a, b))
            bestobjval = zero(T)
            for r in allrays(prob)
                objval = get(lpm.obj) ⋅ r
                if _better(objval, bestobjval)
                    bestobjval = objval
                    lpm.solution = coord(r)
                end
            end
            if _better(bestobjval, zero(T))
                lpm.status = :Unbounded
                lpm.objval = lpm.sense == :Max ? typemax(T) : typemin(T)
            end
            if lpm.status != :Unbounded
                for p in points(prob)
                    objval = get(lpm.obj) ⋅ p
                    if lpm.status == :Undecided || better(objval, get(lpm.objval))
                        lpm.status = :Optimal
                        lpm.objval = objval
                        lpm.solution = p
                    end
                end
            end
        else
            lpm.status = :Optimal
            lpm.objval = zero(T)
            lpm.solution = spzeros(T,N)
        end
    end
end

function MPB.status(lpm::VRepPolyhedraModel)
    lpm.status
end
function MPB.getobjval(lpm::VRepPolyhedraModel)
    get(lpm.objval)
end
function MPB.getsolution(lpm::VRepPolyhedraModel)
    copy(get(lpm.solution))
end
function MPB.getunboundedray(lpm::VRepPolyhedraModel)
    copy(get(lpm.solution))
end
