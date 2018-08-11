# wrapper to convert LPQP solver into Polyhedra solver

# To enable Polyhedra support from an LPQP solver, define, e.g.,
# PolyhedraModel(s::GurobiSolver) = LPQPtoPolyhedraBridge(LinearQuadraticModel(s))

# Cheating a bit here
#PolyhedraModel(s::MathProgbase.AbstractMathProgSolver) = LPQPtoPolyhedraBridge(LinearQuadraticModel(s))
PolyhedraModel(s) = LPQPtoPolyhedraBridge(MPBSI.LinearQuadraticModel(s))

struct LPQPtoPolyhedraBridge <: AbstractPolyhedraModel
    lpqpmodel::MPB.AbstractLinearQuadraticModel
    rep
    c
    sense
end

LPQPtoPolyhedraBridge(m::MPB.AbstractLinearQuadraticModel) = LPQPtoPolyhedraBridge(m, nothing, nothing, nothing)

export LPQPtoPolyhedraBridge

# To transform Polyhedra problems into LinearQuadratic problems
function MPBSI.loadproblem!(m::LPQPtoPolyhedraBridge, rep::HRep, c, sense)
    lp = LPHRepresentation(rep)
    MPBSI.loadproblem!(m.lpqpmodel, lp.A, lp.l, lp.u, c, lp.lb, lp.ub, sense)
end

for f in [:optimize!, :status, :getsolution, :getobjval, :getreducedcosts, :getconstrduals]
    @eval MPBSI.$f(model::LPQPtoPolyhedraBridge) = MPBSI.$f(model.lpqpmodel)
end
