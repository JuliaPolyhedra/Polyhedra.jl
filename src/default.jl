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
default_solver(p::Rep) = nothing
function default_solver(p::Rep, ps::Rep...)
    s = default_solver(p)
    if s === nothing
        default_solver(ps...)
    else
        s
    end
end

"""
    solver(p::Rep, solver::MathProgBase.AbstractMathProgSolver=default_solver(p))

If the V-representation of `p` has been computed, returns `VRepSolver()`, otherwise, returns `solver`.
"""
function solver(p::Rep, solver::MPB.AbstractMathProgSolver=default_solver(p))
    if vrepiscomputed(p)
        VRepSolver()
    else
        solver
    end
end

solver(v::VRepresentation, solver::MPB.AbstractMathProgSolver=default_solver(v)) = VRepSolver()
solver(h::HRepresentation, solver::MPB.AbstractMathProgSolver=default_solver(h)) = solver

_promote_reptype(P1::Type{<:HRep}, ::Type{<:HRep}) = P1
_promote_reptype(P1::Type{<:VRep}, ::Type{<:VRep}) = P1
function promote_reptype(P1::Type{<:Rep}, P2::Type{<:Rep})
    _promote_reptype(P1, P2)
end
promote_reptype(P1::Type{<:Rep}, P2::Type{<:Rep}, P::Type{<:Rep}...) = promote_reptype(promote_reptype(P1, P2), P...)
promote_reptype(P::Type{<:Rep}) = P

# whether Rep has the constructor Rep(::It...; solver=...)
supportssolver(::Type{<:Rep}) = false

function constructpolyhedron(RepT::Type{<:Rep{N, T}}, p::Tuple{Vararg{Rep}}, it::It{N, T}...) where {N, T}
    if supportssolver(RepT)
        solver = default_solver(p...)
        if solver !== nothing && !(solver isa JuMP.UnsetSolver)
            return RepT(it..., solver=solver)
        end
    end
    RepT(it...)::RepT # FIXME without this type annotation even convexhull(::PointsHull{2,Int64,Array{Int64,1}}, ::PointsHull{2,Int64,Array{Int64,1}}) is not type stable, why ?

end

function default_similar(p::Tuple{Vararg{Rep}}, d::FullDim{N}, ::Type{T}, it::It{N, T}...) where {N, T}
    # Some types in p may not support `d` or `T` so we call `similar_type` after `promote_reptype`
    RepT = similar_type(promote_reptype(typeof.(p)...), d, T)
    constructpolyhedron(RepT, p, it...)
end

"""
    similar(p::Tuple{Vararg{Polyhedra.Rep}}, ::Polyhedra.FullDim{N}, ::Type{T}, it::Polyhedra.It{N, T}...)

Creates a representation with a type similar to `p` of a polyhedron of full dimension `N`, element type `T` and initialize it with the iterators `it`.
The type of the result will be chosen closer to the type of `p[1]`.
"""
Base.similar(p::Tuple{Vararg{Rep}}, d::FullDim{N}, ::Type{T}, it::It{N, T}...) where {N, T} = default_similar(p, d, T, it...)
function promote_coefficienttype(p::Tuple{Vararg{Rep}})
    promote_type(MultivariatePolynomials.coefficienttype.(p)...)
end
function Base.similar(p::Tuple{Vararg{Rep{N}}}, d::FullDim{N}, it::It{N}...) where N
    T = promote_coefficienttype(p)
    similar(p, d, T, it...)
end
function Base.similar(p::Tuple{Vararg{Rep{N}}}, it::It{N}...) where N
    similar(p, FullDim{N}(), it...)
end
Base.similar(p::Rep, args...) = similar((p,), args...)
