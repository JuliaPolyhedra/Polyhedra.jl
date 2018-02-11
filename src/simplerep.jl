export SimpleHRepresentation, SimpleVRepresentation

# H-Representation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct SimpleHRepresentation{N, T} <: HRepresentation{N, T}
    # Ax <= b
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    linset::IntSet

    function SimpleHRepresentation{N, T}(A::AbstractMatrix, b::AbstractVector, linset::IntSet=IntSet()) where {N, T}
        if size(A, 1) != length(b)
            error("The length of b must be equal to the number of rows of A")
        end
        if !isempty(linset) && last(linset) > length(b)
            error("The elements of linset should be between 1 and the number of rows of A/length of b")
        end
        if size(A, 2) != N
            error("dimension does not match")
        end
        new{N, T}(A, b, linset)
    end
end

similar_type{N,T}(::Type{<:SimpleHRepresentation}, ::FullDim{N}, ::Type{T}) = SimpleHRepresentation{N,T}
arraytype(p::Union{SimpleHRepresentation{N, T}, Type{SimpleHRepresentation{N, T}}}) where {N, T} = Vector{T}

linset(rep::SimpleHRepresentation) = copy(rep.linset)

function SimpleHRepresentation(A::AbstractMatrix{S}, b::AbstractVector{T}, linset::IntSet=IntSet()) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    SimpleHRepresentation{size(A,2),U}(AbstractMatrix{U}(A), AbstractVector{U}(b), linset)
end
SimpleHRepresentation(A::AbstractMatrix{T}, b::AbstractVector{T}, linset::IntSet=IntSet()) where {T <: Real} = SimpleHRepresentation{size(A,2),T}(A, b, linset)

SimpleHRepresentation(h::HRep{N,T}) where {N,T} = SimpleHRepresentation{N,T}(h)

#function SimpleHRepresentation{N, T}(it::HRepIterator{N, T}) where {N, T}
#    A = Matrix{T}(length(it), N)
#    b = Vector{T}(length(it))
#    linset = IntSet()
#    for (i, h) in enumerate(it)
#        A[i,:] = h.a
#        b[i] = h.β
#        if islin(h)
#            push!(linset, i)
#        end
#    end
#    SimpleHRepresentation{N, T}(A, b, linset)
#end

function SimpleHRepresentation{N, T}(hyperplanes::ElemIt{<:HyperPlane{N, T}}, halfspaces::ElemIt{<:HalfSpace{N, T}}) where {N, T}
    nhyperplane = length(hyperplanes)
    nhrep = nhyperplane + length(halfspaces)
    A = Matrix{T}(nhrep, N)
    b = Vector{T}(nhrep)
    linset = IntSet(1:nhyperplane)
    for (i, h) in enumerate(hyperplanes)
        A[i,:] = h.a
        b[i] = h.β
    end
    for (i, h) in enumerate(halfspaces)
        A[nhyperplane+i,:] = h.a
        b[nhyperplane+i] = h.β
    end
    SimpleHRepresentation{N, T}(A, b, linset)
end

Base.copy(ine::SimpleHRepresentation{N,T}) where {N,T} = SimpleHRepresentation{N,T}(copy(ine.A), copy(ine.b), copy(ine.linset))

Base.length(ine::SimpleHRepresentation) = size(ine.A, 1) # TODO remove

nhyperplanes(ine::SimpleHRepresentation) = length(ine.linset)
starthyperplane(ine::SimpleHRepresentation) = start(ine.linset)
donehyperplane(ine::SimpleHRepresentation, state) = done(ine.linset, state)
function nexthyperplane(ine::SimpleHRepresentation, state)
    (i, nextstate) = next(ine.linset, state)
    (HyperPlane(ine.A[i,:], ine.b[i]), nextstate)
end

function nextz(is::IntSet, i)
    while i in is
        i += 1
    end
    i
end
nhalfspaces(ine::SimpleHRepresentation) = size(ine.A, 1) - nhyperplanes(ine)
starthalfspace(ine::SimpleHRepresentation) = nextz(ine.linset, 1)
donehalfspace(ine::SimpleHRepresentation, state) = state > size(ine.A, 1)
nexthalfspace(ine::SimpleHRepresentation, state) = (HalfSpace(ine.A[state,:], ine.b[state]), nextz(ine.linset, state+1))

function filterintset(J::IntSet, I)
    K = IntSet()
    for (idx, i) in enumerate(I)
        if i in J
            push!(K, idx)
        end
    end
    K
end
Base.getindex(h::SimpleHRepresentation, I::AbstractArray) = SimpleHRepresentation(h.A[I, :], h.b[I], filterintset(h.linset, I))

# V-Representation

mutable struct SimpleVRepresentation{N,T} <: VRepresentation{N,T}
    V::AbstractMatrix{T} # each row is a vertex
    R::AbstractMatrix{T} # each row is a ray
    Vlinset::IntSet
    Rlinset::IntSet

    function SimpleVRepresentation{N, T}(V::AbstractMatrix, R::AbstractMatrix, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) where {N, T}
        if iszero(size(V, 1)) && !iszero(size(R, 1))
            vconcistencyerror()
        end
        if (length(R) > 0 && size(R, 2) != N) || (length(V) > 0 && size(V, 2) != N)
            error("dimension does not match")
        end
        if !isempty(Vlinset) && last(Vlinset) > size(V, 1)
            error("The elements of Vlinset should be between 1 and the number of rows of V")
        end
        if !isempty(Rlinset) && last(Rlinset) > size(R, 1)
            error("The elements of Rlinset should be between 1 and the number of rows of R")
        end
        new{N, T}(V, R, Vlinset, Rlinset)
    end
end

similar_type{N,T}(::Type{<:SimpleVRepresentation}, ::FullDim{N}, ::Type{T}) = SimpleVRepresentation{N,T}
arraytype(p::Union{SimpleVRepresentation{N, T}, Type{SimpleVRepresentation{N, T}}}) where {N, T} = Vector{T}

function linset(rep::SimpleVRepresentation)
    ls = copy(rep.Rlinset)
    nr = nrays(rep)
    for i in rep.Vlinset
        push!(ls, nr + i)
    end
    ls
end

function SimpleVRepresentation(V::AbstractMatrix{S}, R::AbstractMatrix{T}, Vlinset::IntSet=IntSet(), Rlinset::IntSet=IntSet()) where {S <: Real, T <: Real}
    U = promote_type(S, T)
    SimpleVRepresentation{size(V,2),U}(AbstractMatrix{U}(V), AbstractMatrix{U}(R), Vlinset, Rlinset)
end
SimpleVRepresentation(V::AbstractMatrix{T}, linset::IntSet=IntSet()) where {T <: Real} = SimpleVRepresentation{size(V, 2),T}(V, similar(V, 0, size(V, 2)), linset, IntSet())

SimpleVRepresentation(v::VRep{N,T}) where {N,T} = SimpleVRepresentation{N,T}(v)

#function SimpleVRepresentation{N, T}(it::VRepIterator{N, T}) where {N, T}
#    A = Matrix{T}(length(it), N)
#    Rlinset = IntSet()
#    Vlinset = IntSet()
#    points = Int[]
#    rays = Int[]
#    for (i, v) in enumerate(it)
#        A[i,:] = coord(v)
#        if isray(v)
#            push!(rays, i)
#            if islin(v)
#                push!(Rlinset, length(rays))
#            end
#        else
#            push!(points, i)
#            if islin(v)
#                push!(Vlinset, length(points))
#            end
#        end
#    end
#    V = A[points, :]
#    R = A[rays, :]
#    SimpleVRepresentation{N, T}(V, R, Vlinset, Rlinset)
#end

function SimpleVRepresentation{N, T}(sympoints::ElemIt{<:SymPoint{N, T}}, points::ElemIt{<:MyPoint{N, T}}, lines::ElemIt{<:Line{N, T}}, rays::ElemIt{<:Ray{N, T}}) where {N, T}
    nsympoint = length(sympoints)
    npoint = length(points)
    nline = length(lines)
    nray = length(rays)
    V = Matrix{T}(nsympoint + npoint, N)
    R = Matrix{T}(nline + nray, N)
    Vlinset = IntSet()
    Rlinset = IntSet()
    function _fill!(M, linset, offset, ps)
        for (i, p) in enumerate(ps)
            M[offset + i, :] = coord(p)
            if islin(p)
                push!(linset, offset + i)
            end
        end
    end
    _fill!(V, Vlinset, 0, sympoints)
    _fill!(V, Vlinset, nsympoint, points)
    _fill!(R, Rlinset, 0, lines)
    _fill!(R, Rlinset, nline, rays)
    SimpleVRepresentation{N, T}(V, R, Vlinset, Rlinset)
end

Base.copy(ext::SimpleVRepresentation{N,T}) where {N,T} = SimpleVRepresentation{N,T}(copy(ext.V), copy(ext.R), copy(ext.Vlinset), copy(ext.Rlinset))

nsympoints(ext::SimpleVRepresentation) = length(ext.Vlinset)
startsympoint(ext::SimpleVRepresentation) = start(ext.Vlinset)
donesympoint(ext::SimpleVRepresentation, state) = done(ext.Vlinset, state)
function nextsympoint(ext::SimpleVRepresentation, state)
    (i, nextstate) = next(ext.Vlinset, state)
    (SymPoint(ext.V[i,:]), nextstate)
end

npoints(ext::SimpleVRepresentation) = size(ext.V, 1)
startpoint(ext::SimpleVRepresentation) = nextz(ext.Vlinset, 1)
donepoint(ext::SimpleVRepresentation, state) = state > size(ext.V, 1)
nextpoint(ext::SimpleVRepresentation, state) = (ext.V[state,:], nextz(ext.Vlinset, state+1))

nlines(ext::SimpleVRepresentation) = length(ext.Rlinset)
startline(ext::SimpleVRepresentation) = start(ext.Rlinset)
doneline(ext::SimpleVRepresentation, state) = done(ext.Rlinset, state)
function nextline(ext::SimpleVRepresentation, state)
    (i, nextstate) = next(ext.Rlinset, state)
    (Line(ext.R[i,:]), nextstate)
end

nrays(ext::SimpleVRepresentation) = size(ext.R, 1) - nlines(ext)
startray(ext::SimpleVRepresentation) = nextz(ext.Rlinset, 1)
doneray(ext::SimpleVRepresentation, state) = state > size(ext.R, 1)
nextray(ext::SimpleVRepresentation, state) = (Ray(ext.R[state,:]), nextz(ext.Rlinset, state+1))

function Base.getindex(h::SimpleVRepresentation, I::AbstractArray)
    Ir = filter(i -> i <= nrays(h), I)
    Ip = filter(i -> i > nrays(h), I) - nrays(h)
    SimpleVRepresentation(h.V[Ip,:], h.R[Ir,:], filterintset(h.Vlinset, I), filterintset(h.Rlinset, I))
end
