export default_type, default_library, similar_library, library, defaultLPsolverfor

# TODO use vecrep instead of simplerep
default_type(::FullDim{N}, ::Type{T}) where {N, T} = SimplePolyhedron{N, T, SimpleHRepresentation{N, T}, SimpleVRepresentation{N, T}}
default_type(::FullDim{1}, ::Type{T}) where T = Interval{T}

default_library(::FullDim, ::Type{T}) where T = SimplePolyhedraLibrary{T}()
default_library(::FullDim{1}, ::Type{T}) where T = IntervalLibrary{T}()

"""
    similar_library(lib::PolyhedraLibrary, d::FullDim{N}, ::Type{T}) where {N, T}

Returns a library that supports polyhedra of full dimension `T` with coefficient type `T`. If `lib` does not support it, this commonly calls `default_library(d, T)`.
"""
function similar_library end

# Shortcuts
similar_library(p::Polyhedron{N}, d::FullDim, ::Type{T}) where {N, T} = similar_library(library(p), d, T)
similar_library(p::Polyhedron{N}, ::Type{T}) where {N, T} = similar_library(p, FullDim{N}(), T)
similar_library(p::Polyhedron{N, T}, d::FullDim) where {N, T} = similar_library(p, d, T)

"""
    library(p::Polyhedron)

Returns the library used by `p`.
"""
function library end

function getlibrary(args...)
    warn("getlibrary is deprecated, use library instead")
    library(args...)
end

function getlibraryfor(args...)
    warn("getlibraryfor is deprecated, use similar_library instead. Note that the dimension `N` now needs to be given as `FullDim{N}()`.")
    similar_library(p)
end


function defaultLPsolverfor(p::Rep{N,T}, solver=JuMP.UnsetSolver()) where {N,T}
    if vrepiscomputed(p)
        SimpleVRepSolver()
    else
        solver
    end
end

defaultLPsolverfor(::VRepresentation, solver=JuMP.UnsetSolver()) = SimpleVRepSolver()
defaultLPsolverfor(::HRepresentation, solver=JuMP.UnsetSolver()) = solver
