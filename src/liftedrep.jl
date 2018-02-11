export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{N, T} <: HRepresentation{N, T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::AbstractMatrix{T}
    linset::IntSet

    function LiftedHRepresentation{N, T}(A::AbstractMatrix, linset::IntSet=IntSet()) where {N, T}
        if !isempty(linset) && last(linset) > size(A, 1)
            error("The elements of linset should be between 1 and the number of rows of A")
        end
        if size(A, 2) != N+1
            error("dimension does not match")
        end
        new{N, T}(A, linset)
    end
end

similar_type{N,T}(::Type{<:LiftedHRepresentation}, ::FullDim{N}, ::Type{T}) = LiftedHRepresentation{N,T}
arraytype(p::Union{LiftedHRepresentation{N, T}, Type{LiftedHRepresentation{N, T}}}) where {N, T} = Vector{T}

linset(rep::LiftedHRepresentation) = copy(rep.linset)

LiftedHRepresentation(A::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedHRepresentation{size(A,2)-1,T}(A, linset)

LiftedHRepresentation(h::HRepresentation{N,T}) where {N,T} = LiftedHRepresentation{N,T}(h)

#function LiftedHRepresentation{N, T}(it::HRepIterator{N, T}) where {N, T}
#    A = Matrix{T}(length(it), N+1)
#    linset = IntSet()
#    for (i, h) in enumerate(it)
#        A[i,2:end] = -h.a
#        A[i,1] = h.β
#        if islin(h)
#            push!(linset, i)
#        end
#    end
#    LiftedHRepresentation{N, T}(A, linset)
#end

function LiftedHRepresentation{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = Matrix{T}(nhrep, N+1)
    linset = IntSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,2:end] = -h.a
        A[i,1] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,2:end] = -h.a
        A[nhyperplane+i,1] = h.β
    end
    LiftedHRepresentation{N, T}(A, linset)
end

Base.copy(ine::LiftedHRepresentation{N,T}) where {N,T} = LiftedHRepresentation{N,T}(copy(ine.A), copy(ine.linset))

function extractrow(ine::LiftedHRepresentation{N}, i) where N
    β = ine.A[i,1]
    a = -ine.A[i,2:end]
    if i in ine.linset
        HyperPlane(a, β)
    else
        HalfSpace(a, β)
    end
end

nhreps(ine::LiftedHRepresentation) = size(ine.A, 1)

starthrep(ine::LiftedHRepresentation) = 1
donehrep(ine::LiftedHRepresentation, state) = state > nhreps(ine)
nexthrep(ine::LiftedHRepresentation, state) = (extractrow(ine, state), state+1)

nhyperplanes(ine::LiftedHRepresentation) = length(ine.linset)
starthyperplane(ine::LiftedHRepresentation) = start(ine.linset)
donehyperplane(ine::LiftedHRepresentation, state) = done(ine.linset, state)
function nexthyperplane(ine::LiftedHRepresentation{N,T}, state) where {N,T}
    (i, nextstate) = next(ine.linset, state)
    (extractrow(ine, i)::HyperPlane{N,T}, nextstate)
end

nhalfspaces(ine::LiftedHRepresentation) = nhreps(ine) - nhyperplanes(ine)
starthalfspace(ine::LiftedHRepresentation) = nextz(ine.linset, 1)
donehalfspace(ine::LiftedHRepresentation, state) = state > nhreps(ine)
nexthalfspace(ine::LiftedHRepresentation{N,T}, state) where {N,T} = (extractrow(ine, state)::HalfSpace{N,T}, nextz(ine.linset, state+1))

Base.getindex(h::LiftedHRepresentation, I::AbstractArray) = LiftedHRepresentation(h.A[I, :], filterintset(h.linset, I))

# V-Represenation

mutable struct LiftedVRepresentation{N,T} <: VRepresentation{N,T}
    R::AbstractMatrix{T} # each row is a vertex if the first element is 1 and a ray otherwise
    linset::IntSet

    function LiftedVRepresentation{N, T}(R::AbstractMatrix, linset::IntSet=IntSet([])) where {N, T}
        if length(R) > 0 && size(R, 2) != N+1
            error("dimension does not match")
        end
        if !isempty(linset) && last(linset) > size(R, 1)
            error("The elements of linset should be between 1 and the number of rows of R")
        end
        new{N, T}(R, linset)
    end
end

similar_type{N,T}(::Type{<:LiftedVRepresentation}, ::FullDim{N}, ::Type{T}) = LiftedVRepresentation{N,T}
arraytype(p::Union{LiftedVRepresentation{N, T}, Type{LiftedVRepresentation{N, T}}}) where {N, T} = Vector{T}

function linset(rep::LiftedVRepresentation)
    rep.linset
end

LiftedVRepresentation(R::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedVRepresentation{size(R,2)-1,T}(R, linset)

LiftedVRepresentation(v::VRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(v)

#function LiftedVRepresentation{N, T}(it::VRepIterator{N, T}) where {N, T}
#    R = Matrix{T}(length(it), N+1)
#    linset = IntSet()
#    for (i, v) in enumerate(it)
#        R[i,2:end] = coord(v)
#        if isray(v)
#            R[i,1] = zero(T)
#        else
#            R[i,1] = one(T)
#        end
#        if islin(v)
#            push!(linset, i)
#        end
#    end
#    LiftedVRepresentation{N, T}(R, linset)
#end

function LiftedVRepresentation{N, T}(sympoints::ElemIt{<:SymPoint{N, T}}, points::ElemIt{<:MyPoint{N, T}}, lines::ElemIt{<:Line{N, T}}, rays::ElemIt{<:Ray{N, T}}) where {N, T}
    nsympoint = length(sympoints)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    nvrep = nsympoint + npoint + nline + nray
    R = Matrix{T}(nvrep, N+1)
    linset = IntSet()
    function _fill(offset, z, ps)
        for (i, p) in enumerate(ps)
            R[offset + i,2:end] = coord(p)
            R[offset + i,1] = z
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill(0, one(T), sympoints)
    _fill(nsympoint, one(T), points)
    _fill(nsympoint+npoint, zero(T), lines)
    _fill(nsympoint+npoint+nline, zero(T), rays)
    LiftedVRepresentation{N, T}(R, linset)
end

Base.copy(ext::LiftedVRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(copy(ext.R), copy(ext.linset))

_nvreps(ext::LiftedVRepresentation) = size(ext.R, 1)

function isrowpoint(ext::LiftedVRepresentation{N,T}, i) where {N,T}
    ispoint = ext.R[i,1]
    @assert ispoint == zero(T) || ispoint == one(T)
    ispoint == one(T)
end
function _count(ext::LiftedVRepresentation, point::Bool, lin::Bool)
    count = 0
    for i in 1:_nvreps(ext)
        if isrowpoint(ext, i) == point && (i in ext.linset) == lin
            count += 1
        end
    end
    count
end
nsympoints(ext::LiftedVRepresentation) = _count(ext, true, true)
npoints(ext::LiftedVRepresentation) = _count(ext, true, false)
nlines(ext::LiftedVRepresentation) = _count(ext, false, true)
nrays(ext::LiftedVRepresentation) = _count(ext, false, false)

function nextidx(ext::LiftedVRepresentation, i, point::Bool, lin::Bool)
    n = _nvreps(ext)
    while i <= n && (isrowpoint(ext, i) != point || (i in ext.linset) != lin)
        i += 1
    end
    i
end

startsympoint(ext::LiftedVRepresentation) = nextidx(ext, 1, true, true)
donesympoint(ext::LiftedVRepresentation, state) = state > _nvreps(ext)
nextsympoint(ext::LiftedVRepresentation, state) = (SymPoint(ext.R[state,2:end]), nextidx(ext, state+1, true, true))

startpoint(ext::LiftedVRepresentation) = nextidx(ext, 1, true, false)
donepoint(ext::LiftedVRepresentation, state) = state > _nvreps(ext)
nextpoint(ext::LiftedVRepresentation, state) = (ext.R[state,2:end], nextidx(ext, state+1, true, false))

startline(ext::LiftedVRepresentation) = nextidx(ext, 1, false, true)
doneline(ext::LiftedVRepresentation, state) = state > _nvreps(ext)
nextline(ext::LiftedVRepresentation, state) = (Line(ext.R[state,2:end]), nextidx(ext, state+1, false, true))

startray(ext::LiftedVRepresentation) = nextidx(ext, 1, false, false)
doneray(ext::LiftedVRepresentation, state) = state > _nvreps(ext)
nextray(ext::LiftedVRepresentation, state) = (Ray(ext.R[state,2:end]), nextidx(ext, state+1, false, false))

Base.getindex(v::LiftedVRepresentation, I::AbstractArray) = LiftedVRepresentation(v.R[I, :], filterintset(v.linset, I))
