export LiftedHRepresentation, LiftedVRepresentation

# H-Represenation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LiftedHRepresentation{N, T} <: HRepresentation{N, T}
    # Ax >= 0, it is [b -A] * [z; x] where z = 1
    A::AbstractMatrix{T}
    linset::IntSet

    function LiftedHRepresentation{N, T}(A::AbstractMatrix, linset::IntSet=IntSet()) where {N, T}
        if ~isempty(linset) && last(linset) > size(A, 1)
            error("The elements of linset should be between 1 and the number of rows of A")
        end
        if size(A, 2) != N+1
            error("dimension does not match")
        end
        new{N, T}(A, linset)
    end
end

changeeltype{N,T,S}(::Type{LiftedHRepresentation{N,T}}, ::Type{S})  = LiftedHRepresentation{N,S}
changefulldim{N,T}(::Type{LiftedHRepresentation{N,T}}, M)  = LiftedHRepresentation{M,T}
changeboth{N,T,S}(::Type{LiftedHRepresentation{N,T}}, M, ::Type{S}) = LiftedHRepresentation{M,S}

decomposedfast(rep::LiftedHRepresentation) = false
linset(rep::LiftedHRepresentation) = copy(rep.linset)

LiftedHRepresentation(A::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedHRepresentation{size(A,2)-1,T}(A, linset)

LiftedHRepresentation(h::HRepresentation{N,T}) where {N,T} = LiftedHRepresentation{N,T}(h)

function LiftedHRepresentation{N, T}(it::HRepIterator{N, T}) where {N, T}
    A = Matrix{T}(length(it), N+1)
    linset = IntSet()
    for (i, h) in enumerate(it)
        A[i,2:end] = -h.a
        A[i,1] = h.β
        if islin(h)
            push!(linset, i)
        end
    end
    LiftedHRepresentation{N, T}(A, linset)
end

function LiftedHRepresentation{N, T}(eqs, ineqs) where {N, T}
    neq = length(eqs)
    nhrep = neq + length(ineqs)
    A = Matrix{T}(nhrep, N+1)
    linset = IntSet(1:neq)
    for (i, h) in enumerate(eqs)
        A[i,2:end] = -h.a
        A[i,1] = h.β
    end
    for (i, h) in enumerate(ineqs)
        A[neq+i,2:end] = -h.a
        A[neq+i,1] = h.β
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

neqs(ine::LiftedHRepresentation) = length(ine.linset)
starteq(ine::LiftedHRepresentation) = start(ine.linset)
doneeq(ine::LiftedHRepresentation, state) = done(ine.linset, state)
function nexteq(ine::LiftedHRepresentation{N,T}, state) where {N,T}
    (i, nextstate) = next(ine.linset, state)
    (extractrow(ine, i)::HyperPlane{N,T}, nextstate)
end

nineqs(ine::LiftedHRepresentation) = nhreps(ine) - neqs(ine)
startineq(ine::LiftedHRepresentation) = nextz(ine.linset, 1)
doneineq(ine::LiftedHRepresentation, state) = state > nhreps(ine)
nextineq(ine::LiftedHRepresentation{N,T}, state) where {N,T} = (extractrow(ine, state)::HalfSpace{N,T}, nextz(ine.linset, state+1))

Base.getindex(h::LiftedHRepresentation, I::AbstractArray) = LiftedHRepresentation(h.A[I, :], filterintset(h.linset, I))

# V-Represenation

mutable struct LiftedVRepresentation{N,T} <: VRepresentation{N,T}
    R::AbstractMatrix{T} # each row is a vertex if the first element is 1 and a ray otherwise
    linset::IntSet

    function LiftedVRepresentation{N, T}(R::AbstractMatrix, linset::IntSet=IntSet([])) where {N, T}
        if length(R) > 0 && size(R, 2) != N+1
            error("dimension does not match")
        end
        if ~isempty(linset) && last(linset) > size(R, 1)
            error("The elements of linset should be between 1 and the number of rows of R")
        end
        new{N, T}(R, linset)
    end
end

changeeltype{N,T,S}(::Type{LiftedVRepresentation{N,T}}, ::Type{S})  = LiftedVRepresentation{N,S}
changefulldim{N,T}(::Type{LiftedVRepresentation{N,T}}, M)  = LiftedVRepresentation{M,T}
changeboth{N,T,S}(::Type{LiftedVRepresentation{N,T}}, M, ::Type{S}) = LiftedVRepresentation{M,S}

decomposedfast(rep::LiftedVRepresentation) = true
function linset(rep::LiftedVRepresentation)
    rep.linset
end

LiftedVRepresentation(R::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = LiftedVRepresentation{size(R,2)-1,T}(R, linset)

LiftedVRepresentation(v::VRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(v)

function LiftedVRepresentation{N, T}(it::VRepIterator{N, T}) where {N, T}
    R = Matrix{T}(length(it), N+1)
    linset = IntSet()
    for (i, v) in enumerate(it)
        R[i,2:end] = coord(v)
        if isray(v)
            R[i,1] = zero(T)
        else
            R[i,1] = one(T)
        end
        if islin(v)
            push!(linset, i)
        end
    end
    LiftedVRepresentation{N, T}(R, linset)
end

function LiftedVRepresentation{N, T}(points, rays) where {N, T}
    npoint = length(points)
    nray = length(rays)
    nvrep = npoint + nray
    R = Matrix{T}(nvrep, N+1)
    linset = IntSet()
    for (i, p) in enumerate(points)
        R[i,2:end] = coord(p)
        R[i,1] = one(T)
        if islin(p)
            push!(linset, i)
        end
    end
    for (i, r) in enumerate(rays)
        j = npoint + i
        R[j,2:end] = coord(r)
        R[j,1] = zero(T)
        if islin(r)
            push!(linset, j)
        end
    end
    LiftedVRepresentation{N, T}(R, linset)
end

Base.copy(ext::LiftedVRepresentation{N,T}) where {N,T} = LiftedVRepresentation{N,T}(copy(ext.R), copy(ext.linset))

nvreps(ext::LiftedVRepresentation) = size(ext.R, 1)

function isrowpoint(ext::LiftedVRepresentation{N,T}, i) where {N,T}
    ispoint = ext.R[i,1]
    @assert ispoint == zero(T) || ispoint == one(T)
    ispoint == one(T)
end
function extractrow(ext::LiftedVRepresentation{N,T}, i) where {N,T}
    ispoint = ext.R[i,1]
    @assert ispoint == zero(T) || ispoint == one(T)
    a = ext.R[i,2:end]
    if isrowpoint(ext, i)
        if i in ext.linset
            SymPoint(a)
        else
            a
        end
    else
        if i in ext.linset
            Line(a)
        else
            Ray(a)
        end
    end
end
function npoints(ext::LiftedVRepresentation)
    count = 0
    for i in 1:nvreps(ext)
        if isrowpoint(ext, i)
            count += 1
        end
    end
    count
end
nrays(ext::LiftedVRepresentation) = nvreps(ext) - npoints(ext)

startvrep(ext::LiftedVRepresentation) = 1
donevrep(ext::LiftedVRepresentation, state) = state > nvreps(ext)
function nextvrep(ext::LiftedVRepresentation, state)
    (extractrow(ext, state), state+1)
end

function nextrayidx(ext::LiftedVRepresentation, i)
    n = nvreps(ext)
    while i <= n && isrowpoint(ext, i)
        i += 1
    end
    i
end
function nextpointidx(ext::LiftedVRepresentation, i)
    n = nvreps(ext)
    while i <= n && !isrowpoint(ext, i)
        i += 1
    end
    i
end

startray(ext::LiftedVRepresentation) = nextrayidx(ext, 1)
doneray(ext::LiftedVRepresentation, state) = state > nvreps(ext)
function nextray(ext::LiftedVRepresentation, state)
    r = ext.R[state,:]
    (extractrow(ext, state), nextrayidx(ext, state+1))
end

startpoint(ext::LiftedVRepresentation) = nextpointidx(ext, 1)
donepoint(ext::LiftedVRepresentation, state) = state > nvreps(ext)
function nextpoint(ext::LiftedVRepresentation, state)
    p = ext.R[state,:]
    (extractrow(ext, state), nextpointidx(ext, state+1))
end

Base.getindex(v::LiftedVRepresentation, I::AbstractArray) = LiftedVRepresentation(v.R[I, :], filterintset(v.linset, I))
