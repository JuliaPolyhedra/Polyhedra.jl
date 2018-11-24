# wrapper to convert Polyhedra solver into LPQP solver

# To enable LPQP support from a Polyhedra solver, define, e.g.,
# LinearQuadraticModel(s::CDDSolver) = PolyhedraToLPQPBridge(PolyhedraModel(s))

mutable struct PolyhedraToLPQPBridge{T<:Real} <: MPB.AbstractLinearQuadraticModel
    m::AbstractPolyhedraModel
    A::SparseMatrixCSC{T,Int}
    collb::Vector{T}
    colub::Vector{T}
    obj::Vector{T}
    rowlb::Vector{T}
    rowub::Vector{T}
    sense::Symbol
    rowoffset
    coloffset
end

PolyhedraToLPQPBridge(s::AbstractPolyhedraModel) = PolyhedraToLPQPBridge(s, sparse(Int[],Int[],Float64[]), Float64[], Float64[], Float64[], Float64[], Float64[], :Min, nothing, nothing)

export PolyhedraToLPQPBridge

# Loads the provided problem data to set up the linear programming problem:
# min c'x
# st  lb <= Ax <= ub
#      l <=  x <= u
# where sense = :Min or :Max
function MPBSI.loadproblem!(wrap::PolyhedraToLPQPBridge, A, collb, colub, obj, rowlb, rowub, sense)
    wrap.A = A
    wrap.collb = collb
    wrap.colub = colub
    wrap.obj = obj
    wrap.rowlb = rowlb
    wrap.rowub = rowub
    wrap.sense = sense
end

function MPBSI.addquadconstr!(wrap::PolyhedraToLPQPBridge, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
    error("For polyhedra solvers, no quadratic constraint is supported")
end

function computeoffsets(lp::LPHRepresentation)
    coloffset = Vector{Vector{Int}}(undef, size(lp.A, 2))
    rowoffset = Vector{Vector{Int}}(undef, size(lp.A, 1))
    # Assumes that hyperplanes are in first indices
    # and halfspaces in following indices
    curoffset = 0
    for i in 1:size(lp.A, 2)
        if i in lp.coleqs
            curoffset += 1
            coloffset[i] = [curoffset]
        end
    end
    for i in 1:size(lp.A, 1)
        if i in lp.roweqs
            curoffset += 1
            rowoffset[i] = [curoffset]
        end
    end
    for i in 1:size(lp.A, 2)
        if !(i in lp.coleqs)
            offsets = Int[]
            if i in lp.colleqs
                curoffset += 1
                push!(offsets, curoffset)
            end
            if i in lp.colgeqs
                curoffset += 1
                push!(offsets, -curoffset)
            end
            coloffset[i] = offsets
        end
    end
    for i in 1:size(lp.A, 1)
        if !(i in lp.roweqs)
            offsets = Int[]
            if i in lp.rowleqs
                curoffset += 1
                push!(offsets, curoffset)
            end
            if i in lp.rowgeqs
                curoffset += 1
                push!(offsets, -curoffset)
            end
            rowoffset[i] = offsets
        end
    end
    coloffset, rowoffset
end

function MPBSI.optimize!(wrap::PolyhedraToLPQPBridge)
    A = wrap.A
    collb = wrap.collb
    colub = wrap.colub
    obj = wrap.obj
    rowlb = wrap.rowlb
    rowub = wrap.rowub

    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || error("Unequal lengths for row bounds")

    lp = LPHRepresentation(A, collb, colub, rowlb, rowub)
    wrap.coloffset, wrap.rowoffset = computeoffsets(lp)
    MPBSI.loadproblem!(wrap.m, lp, obj, wrap.sense)
    MPBSI.optimize!(wrap.m)
end

MPBSI.getsolution(wrap::PolyhedraToLPQPBridge) = MPBSI.getsolution(wrap.m)
MPBSI.status(wrap::PolyhedraToLPQPBridge) = MPBSI.status(wrap.m)
MPBSI.getobjval(wrap::PolyhedraToLPQPBridge) = MPBSI.getobjval(wrap.m)

function MPBSI.getconstrsolution(wrap::PolyhedraToLPQPBridge)
    m = length(wrap.rowoffset)
    constrsol = MPBSI.getconstrsolution(wrap.m)
    map(offset -> sign(offset[1]) * constrsol[abs(offset[1])], wrap.rowoffset)
end

function getdualsaux(wrap::PolyhedraToLPQPBridge{T}, offsets) where T
    l = length(offsets)
    duals = zeros(T, l)
    constrduals = MPBSI.getconstrduals(wrap.m)
    for i in 1:l
        for offset in offsets[i]
            duals[i] = sign(offset) * constrduals[abs(offset)]
        end
    end
    duals
end

function MPBSI.getreducedcosts(wrap::PolyhedraToLPQPBridge)
    getdualsaux(wrap, wrap.coloffset)
end

function MPBSI.getconstrduals(wrap::PolyhedraToLPQPBridge)
    getdualsaux(wrap, wrap.rowoffset)
end

MPBSI.numconstr(wrap::PolyhedraToLPQPBridge) = size(wrap.A, 1)
MPBSI.numvar(wrap::PolyhedraToLPQPBridge) = size(wrap.A, 2)
MPBSI.getvarLB(wrap::PolyhedraToLPQPBridge) = wrap.collb
MPBSI.getvarUB(wrap::PolyhedraToLPQPBridge) = wrap.colub
MPBSI.getconstrLB(wrap::PolyhedraToLPQPBridge) = wrap.rowlb
MPBSI.getconstrUB(wrap::PolyhedraToLPQPBridge) = wrap.rowub
MPBSI.getobj(wrap::PolyhedraToLPQPBridge) = wrap.obj
MPBSI.getsense(wrap::PolyhedraToLPQPBridge) = wrap.sense
function MPBSI.setvarLB!(wrap::PolyhedraToLPQPBridge, l)
    wrap.collb = l
end
function MPBSI.setvarUB!(wrap::PolyhedraToLPQPBridge, u)
    wrap.colub = u
end
function MPBSI.setconstrLB!(wrap::PolyhedraToLPQPBridge, lb)
    wrap.rowlb = lb
end
function MPBSI.setconstrUB!(wrap::PolyhedraToLPQPBridge, ub)
    wrap.rowub = ub
end
function MPBSI.setobj!(wrap::PolyhedraToLPQPBridge, obj)
    wrap.obj = obj
end
function MPBSI.setsense!(wrap::PolyhedraToLPQPBridge, sense)
    wrap.sense = sense
end
function MPBSI.addvar!(wrap::PolyhedraToLPQPBridge, constridx::AbstractArray{T}, constrcoef, l, u, objcoef) where T<:Integer
    wrap.A = [wrap.A sparsevec(constridx, constrcoef, size(wrap.A, 1))]
    push!(wrap.collb, l)
    push!(wrap.colub, u)
    push!(wrap.obj, objcoef)
end
function MPBSI.addconstr!(wrap::PolyhedraToLPQPBridge, varidx::AbstractArray{T}, coef, lb, ub) where T<:Integer
    wrap.A = [wrap.A; sparsevec(varidx, coef, size(wrap.A, 2))']
    push!(wrap.rowlb, lb)
    push!(wrap.rowub, ub)
end
