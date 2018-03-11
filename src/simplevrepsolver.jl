export SimpleVRepSolver, SimpleVRepPolyhedraModel

struct SimpleVRepSolver <: MathProgBase.AbstractMathProgSolver
end

mutable struct SimpleVRepPolyhedraModel <: AbstractPolyhedraModel
    vrep::Nullable{VRep}
    obj::Nullable{Vector}
    sense::Nullable{Symbol}

    objval::Nullable
    solution::Nullable{Vector}
    status::Symbol

    function SimpleVRepPolyhedraModel()
        new(nothing, nothing, nothing, nothing, nothing, :Undecided)
    end
end

PolyhedraModel(::SimpleVRepSolver) = SimpleVRepPolyhedraModel()
LinearQuadraticModel(s::SimpleVRepSolver) = PolyhedraToLPQPBridge(PolyhedraModel(s))

function loadproblem!(lpm::SimpleVRepPolyhedraModel, vrep::VRep, obj, sense)
    if !(sense in [:Max, :Min])
        error("sense should be :Max or :Min")
    end
    lpm.vrep = vrep
    if !isapproxzero(obj)
        lpm.obj = copy(obj)
        lpm.sense = sense
    end
end

function optimize!(lpm::SimpleVRepPolyhedraModel)
    if isnull(lpm.vrep)
        error("No problem loaded")
    end
    prob = get(lpm.vrep)
    N = fulldim(prob)
    T = MultivariatePolynomials.coefficienttype(prob)
    if !hassympoints(prob) && !haspoints(prob) && !haslines(prob) && !hasrays(prob)
        lpm.status = :Infeasible
    else
        if lpm.sense in [:Max, :Min]
            better(a, b) = (lpm.sense == :Max ? a > b : a < b)
            mybetter(a, b) = (lpm.sense == :Max ? mygt(a, b) : mylt(a, b))
            lpm.status = :Infeasible
            for r in allrays(prob)
                objval = lpm.obj ⋅ r
                if lpm.status != :Unbounded && mybetter(objval, zero(T))
                    lpm.status = :Unbounded
                    lpm.objval = lpm.sense == :Max ? typemax(T) : typemin(T)
                    lpm.solution = r
                end
            end
            if status != :Unbounded
                for p in allpoints(prob)
                    objval = lpm.obj ⋅ p
                    if lpm.status == :Undecided || (lpm.status == :Optimal && better(objval, get(lpm.objval)))
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

function MathProgBase.status(lpm::SimpleVRepPolyhedraModel)
    lpm.status
end
function MathProgBase.getobjval(lpm::SimpleVRepPolyhedraModel)
    get(lpm.objval)
end
function MathProgBase.getsolution(lpm::SimpleVRepPolyhedraModel)
    copy(get(lpm.solution))
end
function MathProgBase.getunboundedray(lpm::SimpleVRepPolyhedraModel)
    copy(get(lpm.solution))
end
