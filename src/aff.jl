export HAffineSpace, VAffineSpace, affinehull, linespace, detecthlinearities!, detectvlinearities!

# Linearity
detectvlinearities!(p::VRep) = error("detectvlinearities! not implemented for $(typeof(p))")
detecthlinearities!(p::HRep) = error("detecthlinearities! not implemented for $(typeof(p))")
function dim{N}(h::HRep{N}, current=false)
    if !current
        detecthlinearities!(h::HRep)
    end
    N - neqs(h)
end


# It is easy to go from H-rep of affine space to V-rep of affine space by computing the kernel of a matrix using RowEchelon
# However, it is really worth it since Base.in for an HRepElement in HAffineSpace and Base.in for an VRepElement in VAffineSpace are false already.
function remproj(x, l)
    x * dot(coord(l), coord(l)) - l * dot(coord(x), coord(l))
end

# H-representation

# Representation of an affine space as the intersection of hyperplanes.
# Also called affine set, affine manifold, affine variety, linear variety or flat.
# An affine space L satisfies:
# λx + (1-λ)y ∈ L, ∀x, y ∈ L, ∀ λ ∈ R
# Note that λ is not required to be between 0 and 1 as in convex sets.
immutable HAffineSpace{N, T} <: HRepresentation{N, T}
    # HyperPlanes whose intersection is the affine space
    hps::Vector{HyperPlane{N, T}}
end

(::Type{HAffineSpace{N, T}}){N, T}() = HAffineSpace{N, T}(HyperPlane{N, T}[])
function (::Type{HAffineSpace{N, T}}){N, T}(it::EqIterator)
    HAffineSpace{N, T}([hp for hp in it])
end
HAffineSpace{N, T}(it::EqIterator{N, T}) = HAffineSpace{N, T}(it)

Base.push!{N, T}(L::HAffineSpace{N, T}, h::HyperPlane{N, T}) = push!(L.hps, h)

decomposedfast(L::HAffineSpace) = true

nineqs(::HAffineSpace) = 0
startineq(::HAffineSpace) = 0
doneineq(::HAffineSpace, state) = true
neqs(L::HAffineSpace) = length(L.hps)
starteq(L::HAffineSpace) = start(L.hps)
nexteq(L::HAffineSpace, state) = next(L.hps, state)
doneeq(L::HAffineSpace, state) = done(L.hps, state)

# Returns an HAffineSpace representing the affine hull of p.
# The affine hull is defined as
# {λx + (1-λ)y | x, y ∈ p, λ ∈ R}
function affinehull(h::HRep, current=false)
    if !current
        detecthlinearities!(h)
    end
    HAffineSpace(eqs(h))
end

function remproj(h::HRepElement, L::HAffineSpace)
    for hp in eqs(L)
        h = remproj(h, hp)
    end
    h
end
function Base.in(h::HRepElement, L::HAffineSpace)
    h = remproj(h, L)
    myeqzero(h)
end

function removeduplicates{N, T}(L::HAffineSpace{N, T})
    H = HAffineSpace{N, T}()
    for h in eqs(L)
        if !(h in H)
            push!(H, h)
        end
    end
    H
end

# V-representation

# Representation of an affine space containing the origin by the minkowsky sum of lines
immutable VAffineSpace{N, T} <: VRepresentation{N, T}
    lines::Vector{Line{N, T}}
end

(::Type{VAffineSpace{N, T}}){N, T}() = VAffineSpace{N, T}(Line{N, T}[])
function (::Type{VAffineSpace{N, T}}){N, T}(it::LineIterator)
    VAffineSpace{N, T}([l for l in it])
end
VAffineSpace{N, T}(it::LineIterator{N, T}) = VAffineSpace{N, T}(it)

Base.push!{N, T}(L::VAffineSpace{N, T}, l::Line{N, T}) = push!(L.lines, l)

decomposedfast(L::VAffineSpace) = true

npoint(L::VAffineSpace) = 0
startpoint(L::VAffineSpace) = 0
donepoint(L::VAffineSpace, state) = true
nrays(L::VAffineSpace) = nlines(L)
startray(L::VAffineSpace) = start(L.lines)
doneray(L::VAffineSpace, state) = done(L.lines, state)
nextray(L::VAffineSpace, state) = next(L.lines, state)
nlines(L::VAffineSpace) = length(L.lines)
startline(L::VAffineSpace) = start(L.lines)
doneline(L::VAffineSpace, state) = done(L.lines, state)
nextline(L::VAffineSpace, state) = next(L.lines, state)

# Returns a VAffineSpace representing the following set (TODO does it have a name?, does someone has a reference talking about it ?)
# {x | ⟨a, x⟩ = 0 ∀ a such that (α, β) is a valid inequality for p}
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

function removeduplicates{N, T}(L::VAffineSpace{N, T})
    V = VAffineSpace{N, T}()
    for h in eqs(L)
        if !(h in H)
            push!(H, h)
        end
    end
    H
end
