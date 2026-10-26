export DefaultLibrary, DefaultPolyhedron

"""
    DefaultLibrary{T}

Default library for polyhedra of dimension larger than 1 ([`IntervalLibrary`](@ref) is the default for polyhedra of dimension 1).
The library implements the bare minimum and uses the fallback implementation for all operations.
"""
struct DefaultLibrary{T} <: Library
    solver::Any
    function DefaultLibrary{T}(solver = nothing) where {T}
        return new{T}(solver)
    end
end

function similar_library(lib::DefaultLibrary, d::FullDim, ::Type{T}) where {T}
    return default_library(d, T)
end # default_library allows to fallback to Interval if d is FullDim{1}

mutable struct DefaultPolyhedron{
    T,
    HRepT<:HRepresentation{T},
    VRepT<:VRepresentation{T},
} <: Polyhedron{T}
    hrep::Union{HRepT,Nothing}
    hred::Redundancy
    vrep::Union{VRepT,Nothing}
    vred::Redundancy
    solver::Any
    function DefaultPolyhedron{T,HRepT,VRepT}(
        hrep::Union{HRepT,Nothing},
        vrep::Union{VRepT,Nothing},
        solver,
    ) where {T,HRepT<:HRepresentation{T},VRepT<:VRepresentation{T}}
        return new{T,HRepT,VRepT}(
            hrep,
            UNKNOWN_REDUNDANCY,
            vrep,
            UNKNOWN_REDUNDANCY,
            solver,
        )
    end
end
function DefaultPolyhedron{T,HRepT,VRepT}(
    hrep::HRepT,
    solver,
) where {T,HRepT<:HRepresentation{T},VRepT<:VRepresentation{T}}
    return DefaultPolyhedron{T,HRepT,VRepT}(hrep, nothing, solver)
end
function DefaultPolyhedron{T,HRepT,VRepT}(
    vrep::VRepT,
    solver,
) where {T,HRepT<:HRepresentation{T},VRepT<:VRepresentation{T}}
    return DefaultPolyhedron{T,HRepT,VRepT}(nothing, vrep, solver)
end
function DefaultPolyhedron{T,HRepT,VRepT}(
    hrep::HRepresentation,
    solver,
) where {T,HRepT<:HRepresentation{T},VRepT<:VRepresentation{T}}
    return DefaultPolyhedron{T,HRepT,VRepT}(convert(HRepT, hrep), solver)
end
function DefaultPolyhedron{T,HRepT,VRepT}(
    vrep::VRepresentation,
    solver,
) where {T,HRepT<:HRepresentation{T},VRepT<:VRepresentation{T}}
    return DefaultPolyhedron{T,HRepT,VRepT}(convert(VRepT, vrep), solver)
end

FullDim(p::DefaultPolyhedron) = FullDim_rep(p.hrep, p.vrep)
library(::Type{<:DefaultPolyhedron{T}}) where {T} = DefaultLibrary{T}()
library(p::DefaultPolyhedron{T}) where {T} = DefaultLibrary{T}(p.solver)
default_solver(p::DefaultPolyhedron; T = nothing) = p.solver
supportssolver(::Type{<:DefaultPolyhedron}) = true

function hvectortype(::Type{<:DefaultPolyhedron{T,HRepT}}) where {T,HRepT}
    return hvectortype(HRepT)
end
function vvectortype(
    ::Type{DefaultPolyhedron{T,HRepT,VRepT}},
) where {T,HRepT,VRepT}
    return vvectortype(VRepT)
end

function similar_type(
    ::Type{<:DefaultPolyhedron{S,HRepT,VRepT}},
    d::FullDim,
    ::Type{T},
) where {S,HRepT,VRepT,T}
    return DefaultPolyhedron{
        T,
        similar_type(HRepT, d, T),
        similar_type(VRepT, d, T),
    }
end

function DefaultPolyhedron{T,HRepT,VRepT}(
    d::FullDim,
    hits::HIt...;
    solver = nothing,
) where {T,HRepT,VRepT}
    return DefaultPolyhedron{T,HRepT,VRepT}(HRepT(d, hits...), solver)
end
function DefaultPolyhedron{T,HRepT,VRepT}(
    d::FullDim,
    vits::VIt...;
    solver = nothing,
) where {T,HRepT,VRepT}
    return DefaultPolyhedron{T,HRepT,VRepT}(VRepT(d, vits...), solver)
end

# Need fulltype in case the use does `intersect!` with another element
function DefaultPolyhedron{T}(rep::Representation, solver) where {T}
    return DefaultPolyhedron{T}(change_coefficient_type(rep, T), solver)
end
function DefaultPolyhedron{T}(rep::HRepresentation{T}, solver) where {T}
    HRepT = fulltype(typeof(rep))
    VRepT = dualtype(HRepT)
    return DefaultPolyhedron{T,HRepT,VRepT}(rep, solver)
end
function DefaultPolyhedron{T}(rep::VRepresentation{T}, solver) where {T}
    VRepT = fulltype(typeof(rep))
    HRepT = dualtype(VRepT)
    return DefaultPolyhedron{T,HRepT,VRepT}(rep, solver)
end

function polyhedron(rep::Representation, lib::DefaultLibrary{T}) where {T}
    return DefaultPolyhedron{polytypefor(T)}(rep, lib.solver)
end

function Base.copy(p::DefaultPolyhedron{T}) where {T}
    if p.hrep !== nothing
        DefaultPolyhedron{T}(p.hrep, p.solver)
    else
        DefaultPolyhedron{T}(p.vrep, p.solver)
    end
end

hredundancy(p::DefaultPolyhedron) = p.hred
vredundancy(p::DefaultPolyhedron) = p.vred

hrepiscomputed(p::DefaultPolyhedron) = p.hrep !== nothing
function computehrep!(p::DefaultPolyhedron)
    # vrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.hrep = doubledescription(p.vrep, p.solver)
    return p.hred = NO_REDUNDANCY
end
function hrep(p::DefaultPolyhedron)
    if !hrepiscomputed(p)
        computehrep!(p)
    end
    return p.hrep
end
vrepiscomputed(p::DefaultPolyhedron) = p.vrep !== nothing
function computevrep!(p::DefaultPolyhedron)
    # hrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.vrep = doubledescription(p.hrep, p.solver)
    return p.vred = NO_REDUNDANCY
end
function vrep(p::DefaultPolyhedron)
    if !vrepiscomputed(p)
        computevrep!(p)
    end
    return p.vrep
end

# Fallback to avoid breaking library adding this new `red` arg.
# TODO We should remove these fallbacks in the future
sethrep!(p, h, red) = sethrep!(p, h)
setvrep!(p, v, red) = setvrep!(p, v)
resethrep!(p, h, red) = resethrep!(p, h)
resetvrep!(p, v, red) = resetvrep!(p, v)

function sethrep!(
    p::DefaultPolyhedron,
    h::HRepresentation,
    redundancy::Redundancy = UNKNOWN_REDUNDANCY,
)
    p.hrep = h
    p.hred = redundancy
    return
end
function setvrep!(
    p::DefaultPolyhedron,
    v::VRepresentation,
    redundancy::Redundancy = UNKNOWN_REDUNDANCY,
)
    p.vrep = v
    p.vred = redundancy
    return
end
function resethrep!(
    p::DefaultPolyhedron,
    h::HRepresentation,
    redundancy::Redundancy = UNKNOWN_REDUNDANCY,
)
    p.hrep = h
    p.hred = redundancy
    p.vrep = nothing
    return
end
function resetvrep!(
    p::DefaultPolyhedron,
    v::VRepresentation,
    redundancy::Redundancy = UNKNOWN_REDUNDANCY,
)
    p.vrep = v
    p.vred = redundancy
    p.hrep = nothing
    return
end
