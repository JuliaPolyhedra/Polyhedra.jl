export HAffineSpace, VAffineSpace, affinehull, linespace, detecthlinearities!, detectvlinearities!

# Linearity
detectvlinearities!(p::VRep) = error("detectvlinearities! not implemented for $(typeof(p))")
detecthlinearities!(p::HRep) = error("detecthlinearities! not implemented for $(typeof(p))")
function dim(h::HRep{N}, current=false) where N
    if !current
        detecthlinearities!(h::HRep)
    end
    N - nhyperplanes(h)
end


# It is easy to go from H-rep of affine space to V-rep of affine space by computing the kernel of a matrix using RowEchelon
# However, it is really worth it since Base.in for an HRepElement in HAffineSpace and Base.in for an VRepElement in VAffineSpace are false already.
function remproj(x::RepElement{N, <:Integer}, l::RepElement{N, <:Integer}) where N
    # generates large numbers but keeps the type integer
    x * dot(coord(l), coord(l)) - l * dot(coord(x), coord(l))
end
function remproj(x, l)
    simplify(x - l * (dot(coord(x), coord(l)) / dot(coord(l), coord(l))))
end

# H-representation

# Representation of an affine space as the intersection of hyperplanes.
# Also called affine set, affine manifold, affine variety, linear variety or flat.
# An affine space L satisfies:
# λx + (1-λ)y ∈ L, ∀x, y ∈ L, ∀ λ ∈ R
# Note that λ is not rhyperplaneuired to be between 0 and 1 as in convex sets.
struct HAffineSpace{N, T, AT} <: HRepresentation{N, T}
    # HyperPlanes whose intersection is the affine space
    hps::Vector{HyperPlane{N, T, AT}}
    function HAffineSpace{N, T, AT}(hps::Vector{HyperPlane{N, T, AT}}) where {N, T, AT}
        new{N, T, AT}(hps)
    end
end
arraytype(L::HAffineSpace{N, T, AT}) where {N, T, AT} = AT

HAffineSpace{N, T, AT}() where {N, T, AT} = HAffineSpace{N, T, AT}(HyperPlane{N, T, AT}[])
function HAffineSpace{N, T, AT}(it::ElemIt{HyperPlane{N, T, AT}}) where {N, T, AT}
    HAffineSpace{N, T, AT}(collect(it))
end
HAffineSpace(it::ElemIt{HyperPlane{N, T, AT}}) where {N, T, AT} = HAffineSpace{N, T, AT}(it)

Base.intersect!(L::HAffineSpace{N}, h::HyperPlane{N}) where N = push!(L.hps, h)

nhalfspaces(::HAffineSpace) = 0
starthalfspace(::HAffineSpace) = 0
donehalfspace(::HAffineSpace, state) = true
nhyperplanes(L::HAffineSpace) = length(L.hps)
starthyperplane(L::HAffineSpace) = start(L.hps)
nexthyperplane(L::HAffineSpace, state) = next(L.hps, state)
donehyperplane(L::HAffineSpace, state) = done(L.hps, state)

# Returns an HAffineSpace representing the affine hull of p.
# The affine hull is defined as
# {λx + (1-λ)y | x, y ∈ p, λ ∈ R}
function affinehull(h::HRep, current=false)
    if !current
        detecthlinearities!(h)
    end
    HAffineSpace(hyperplanes(h))
end

function remproj(h::HRepElement, L::HAffineSpace)
    for hp in hyperplanes(L)
        h = remproj(h, hp)
    end
    h
end
function Base.in(h::HRepElement, L::HAffineSpace)
    h = remproj(h, L)
    myeqzero(h)
end

function removeduplicates(L::HAffineSpace{N, T, AT}) where {N, T, AT}
    H = HAffineSpace{N, T, AT}()
    for h in hyperplanes(L)
        if !(h in H)
            intersect!(H, h)
        end
    end
    H
end

# V-representation

abstract type VCone{N, T, AT} <: VRepresentation{N, T} end

nsympoint(L::VCone) = 0
startsympoint(L::VCone) = true
donesympoint(L::VCone, state::Bool) = state

# See issue #28
npoint(L::VCone) = 1
startpoint(L::VCone) = false
donepoint(L::VCone, state::Bool) = state
function nextpoint(L::VCone{N, T, AT}, state::Bool) where {N, T, AT}
    @assert !state
    zeros(AT)
end

# Representation of an affine space containing the origin by the minkowsky sum of lines
struct VAffineSpace{N, T, AT} <: VCone{N, T, AT}
    lines::Vector{Line{N, T, AT}}
    function VAffineSpace{N, T, AT}(lines::Vector{Line{N, T, AT}}) where {N, T, AT}
        new{N, T, AT}(lines)
    end
end
arraytype(L::VAffineSpace{N, T, AT}) where {N, T, AT} = AT

VAffineSpace{N, T, AT}() where {N, T, AT} = VAffineSpace{N, T, AT}(Line{N, T, AT}[])
function VAffineSpace{N, T, AT}(it::ElemIt{Line{N, T, AT}}) where {N, T, AT}
    VAffineSpace{N, T, AT}(collect(it))
end
VAffineSpace(it::ElemIt{Line{N, T, AT}}) where {N, T, AT} = VAffineSpace{N, T, AT}(it)

convexhull!(L::VAffineSpace{N}, l::Line{N}) where {N} = push!(L.lines, l)

nrays(L::VAffineSpace) = nlines(L)
startray(L::VAffineSpace) = start(L.lines)
doneray(L::VAffineSpace, state) = done(L.lines, state)
nextray(L::VAffineSpace, state) = next(L.lines, state)
nlines(L::VAffineSpace) = length(L.lines)
startline(L::VAffineSpace) = start(L.lines)
doneline(L::VAffineSpace, state) = done(L.lines, state)
nextline(L::VAffineSpace, state) = next(L.lines, state)

# Returns a VAffineSpace representing the following set (TODO does it have a name?, does someone has a reference talking about it ?)
# {x | ⟨a, x⟩ = 0 ∀ a such that (α, β) is a valid inhyperplaneuality for p}
function linespace(v::VRep, current=false)
    if !current
        detectvlinearities!(v)
    end
    VAffineSpace(lines(v))
end

function remproj(v::VRepElement, L::VAffineSpace)
    for l in lines(L)
        v = remproj(v, l)
    end
    v
end
function Base.in(v::VRepElement, L::VAffineSpace)
    v = remproj(v, L)
    myeqzero(coord(v))
end

function removeduplicates(L::VAffineSpace{N, T, AT}) where {N, T, AT}
    V = VAffineSpace{N, T, AT}()
    for h in hyperplanes(L)
        if !(h in H)
            convexhull!(H, h)
        end
    end
    H
end
