# wrapper to convert Polyhedra solver into LPQP solver

# To enable LPQP support from a Polyhedra solver, define, e.g.,
# LinearQuadraticModel(s::CDDSolver) = PolyhedraToLPQPBridge(PolyhedraModel(s))

type PolyhedraToLPQPBridge <: MathProgBase.AbstractLinearQuadraticModel
    m::AbstractPolyhedraModel
    A::SparseMatrixCSC{Float64,Int}
    collb::Vector{Float64}
    colub::Vector{Float64}
    obj::Vector{Float64}
    rowlb::Vector{Float64}
    rowub::Vector{Float64}
    sense::Symbol
end

PolyhedraToLPQPBridge(s::AbstractPolyhedraModel) = PolyhedraToLPQPBridge(s, sparse(Int[],Int[],Float64[]), Float64[], Float64[], Float64[], Float64[], Float64[], :Uninitialized, Int[], Array(Vector{Int},0))

export PolyhedraToLPQPBridge

# Loads the provided problem data to set up the linear programming problem:
# min c'x
# st  lb <= Ax <= ub
#      l <=  x <= u
# where sense = :Min or :Max
function loadproblem!(wrap::PolyhedraToLPQPBridge, A, collb, colub, obj, rowlb, rowub, sense)
    wrap.A = A
    wrap.collb = collb
    wrap.colub = colub
    wrap.obj = obj
    wrap.rowlb = rowlb
    wrap.rowub = rowub
    wrap.sense = sense
end

function addquadconstr!(wrap::PolyhedraToLPQPBridge, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    error("For polyhedra solvers, no quadratic constraint is supported")
end

function optimize!(wrap::PolyhedraToLPQPBridge)
    A = wrap.A
    collb = wrap.collb
    colub = wrap.colub
    obj = wrap.obj
    rowlb = wrap.rowlb
    rowub = wrap.rowub

    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || error("Unequal lengths for row bounds")

    lp = LPHRepresentation(A, collb, colub, rowlb, rowub)
    loadproblem!(wrap.m, obj, A, b, constr_cones, var_cones, wrap.obj, wrap.sense)
    optimize!(wrap.m)
end

getsolution(wrap::PolyhedraToLPQPBridge) = getsolution(wrap.m)
status(wrap::PolyhedraToLPQPBridge) = status(wrap.m)
getobjval(wrap::PolyhedraToLPQPBridge) = getobjval(wrap.m)

function getreducedcosts(wrap::PolyhedraToLPQPBridge)
    warning("Not correct for range constraints")
    getreducedcosts(wrap.m)
end

function getconstrduals(wrap::PolyhedraToLPQPBridge)
    warning("Not correct for range constraints")
    getconstrduals(wrap.m)
end
