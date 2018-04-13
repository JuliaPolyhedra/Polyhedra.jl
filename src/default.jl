export default_type, default_library, similar_library, library
export getlibrary, getlibraryfor

"""
    default_type(::FullDim{N}, ::Type{T}) where {N, T}

Returns the default polyhedron type for `N`-dimensional polyhedron of coefficient type `T`.
"""
function default_type end

default_type(::FullDim{N}, ::Type{T}) where {N, T} = SimplePolyhedron{N, T, Intersection{N, T, Vector{T}}, Hull{N, T, Vector{T}}}
default_type(::FullDim{1}, ::Type{T}) where T = Interval{T, SVector{1, T}}

"""
    default_library(::FullDim{N}, ::Type{T}) where {N, T}

Returns the default polyhedral library for `N`-dimensional polyhedron of coefficient type `T`.
"""
function default_library end

_default_type(::Type{T}) where T = T
# See https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/35
_default_type(::Type{AbstractFloat}) = Float64
default_library(::FullDim, ::Type{T}) where T = SimplePolyhedraLibrary{_default_type(T)}()
default_library(::FullDim{1}, ::Type{T}) where T = IntervalLibrary{_default_type(T)}()

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
    Base.depwarn("getlibrary is deprecated, use library instead", :getlibrary)
    library(args...)
end

function getlibraryfor(args...)
    Base.depwarn("getlibraryfor is deprecated, use similar_library instead. Note that the dimension `N` now needs to be given as `FullDim{N}()`.", :getlibraryfor)
    similar_library(args...)
end


"""
    default_solver(p::Rep)

Returns a default linear programming solver for the polyhedron `p` (e.g. CDD has an internal solver which is used by default).
"""
default_solver(p::Rep) = JuMP.UnsetSolver()

"""
    solver(p::Rep, solver::MathProgBase.AbstractMathProgSolver=default_solver(p))

If the V-representation of `p` has been computed, returns `VRepSolver()`, otherwise, returns `solver`.
"""
function solver(p::Rep, solver::MathProgBase.AbstractMathProgSolver=default_solver(p))
    if vrepiscomputed(p)
        VRepSolver()
    else
        solver
    end
end

solver(v::VRepresentation, solver::MathProgBase.AbstractMathProgSolver=default_solver(v)) = VRepSolver()
solver(h::HRepresentation, solver::MathProgBase.AbstractMathProgSolver=default_solver(h)) = solver
