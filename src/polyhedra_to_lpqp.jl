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

function computeoffsets(lp::LPHRepresentation)
    coloffset = Vector{Vector{Int}}(size(lp.A, 2))
    curoffset = 0
    for i in 1:size(lp.A, 2)
        if i in lp.coleqs
            curoffset += 1
            offsets = [curoffset]
        else
            offsets = Int[]
            if i in lp.colleqs
                curoffset += 1
                push!(offsets, curoffset)
            end
            if i in lp.colgeqs
                curoffset += 1
                push!(offsets, -curoffset)
            end
        end
        coloffset[i] = offsets
    end
    rowoffset = Vector{Vector{Int}}(size(lp.A, 1))
    for i in 1:size(lp.A, 1)
        if i in lp.roweqs
            curoffset += 1
            offsets = [curoffset]
        else
            offsets = Int[]
            if i in lp.rowleqs
                curoffset += 1
                push!(offsets, curoffset)
            end
            if i in lp.rowgeqs
                curoffset += 1
                push!(offsets, -curoffset)
            end
        end
        rowoffset[i] = offsets
    end
    coloffset, rowoffset
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
    wrap.coloffset, wrap.rowoffset = computeoffsets(lp)
    loadproblem!(wrap.m, lp, obj, wrap.sense)
    optimize!(wrap.m)
end

getsolution(wrap::PolyhedraToLPQPBridge) = getsolution(wrap.m)
status(wrap::PolyhedraToLPQPBridge) = status(wrap.m)
getobjval(wrap::PolyhedraToLPQPBridge) = getobjval(wrap.m)

function getconstrsolution(wrap::PolyhedraToLPQPBridge)
    m = length(wrap.rowoffset)
    constrsol = getconstrsolution(wrap.m)
    map(offset -> sign(offset[1]) * constrsol[abs(offset[1])], wrap.rowoffset)
end

function getdualsaux(wrap::PolyhedraToLPQPBridge, offsets)
    l = length(offsets)
    duals = zeros(Float64, l)
    constrduals = getconstrduals(wrap.m)
    for i in 1:l
        for offset in offsets[i]
            duals[i] = sign(offset) * constrduals[abs(offset)]
        end
    end
    duals
end

function getreducedcosts(wrap::PolyhedraToLPQPBridge)
    getdualsaux(wrap, wrap.coloffset)
end

function getconstrduals(wrap::PolyhedraToLPQPBridge)
    getdualsaux(wrap, wrap.rowoffset)
end

numconstr(wrap::PolyhedraToLPQPBridge) = size(wrap.A, 1)
numvar(wrap::PolyhedraToLPQPBridge) = size(wrap.A, 2)
getvarLB(wrap::PolyhedraToLPQPBridge) = wrap.collb
getvarUB(wrap::PolyhedraToLPQPBridge) = wrap.colub
getconstrLB(wrap::PolyhedraToLPQPBridge) = wrap.rowlb
getconstrUB(wrap::PolyhedraToLPQPBridge) = wrap.rowub
getobj(wrap::PolyhedraToLPQPBridge) = wrap.obj
getsense(wrap::PolyhedraToLPQPBridge) = wrap.sense
function setvarLB!(wrap::PolyhedraToLPQPBridge, l)
    wrap.collb = l
end
function setvarUB!(wrap::PolyhedraToLPQPBridge, u)
    wrap.colub = u
end
function setconstrLB!(wrap::PolyhedraToLPQPBridge, lb)
    wrap.rowlb = lb
end
function setconstrUB!(wrap::PolyhedraToLPQPBridge, ub)
    wrap.rowub = ub
end
function setobj!(wrap::PolyhedraToLPQPBridge, obj)
    wrap.obj = obj
end
function setsense!(wrap::PolyhedraToLPQPBridge, sense)
    wrap.sense = sense
end
function addvar!{T<:Integer}(wrap::PolyhedraToLPQPBridge, constridx::AbstractArray{T}, constrcoef, l, u, objcoef)
    wrap.A = [wrap.A sparsevec(constridx, constrcoef, size(wrap.A, 1))]
    push!(wrap.collb, l)
    push!(wrap.colub, u)
    push!(wrap.obj, objcoef)
end
function addconstr!{T<:Integer}(wrap::PolyhedraToLPQPBridge, varidx::AbstractArray{T}, coef, lb, ub)
    wrap.A = [wrap.A; sparsevec(varidx, coef, size(wrap.A, 2))']
    push!(wrap.rowlb, lb)
    push!(wrap.rowub, ub)
end
