export default_type, default_library, similar_library, library

"""
    default_type(d::FullDim, ::Type{T}) where {T}

Returns the default polyhedron type for `d`-dimensional polyhedron of coefficient type `T`.
"""
function default_type end

function default_type(d::StaticArrays.Size{N}, ::Type{T}) where {N, T}
    return DefaultPolyhedron{T, Intersection{T, StaticArrays.SVector{N[1], T}, typeof(d)}, Hull{T, StaticArrays.SVector{N[1], T}, typeof(d)}}
end
function default_type(d::StaticArrays.Size{(1,)}, ::Type{T}) where T
    return Interval{T, StaticArrays.SVector{1, T}, typeof(d)}
end
function default_type(d::Int, ::Type{T}) where T
    if d == 1
        return Interval{T, Vector{T}, typeof(d)}
    else
        return DefaultPolyhedron{T, Intersection{T, Vector{T}, typeof(d)}, Hull{T, Vector{T}, typeof(d)}}
    end
end

"""
    default_library(d::FullDim, ::Type{T}) where {T}

Returns the default polyhedral library for `d`-dimensional polyhedron of coefficient type `T`.

## Examples

To obtain the default library for 2-dimensional polyhedra of eltype `Float64`,
do `default_library(2, Float64)`.

Given an `StaticArrays.SVector` `v`, to obtain a default library for points
of the type of `v` in a type stable way, do
`default_library(Polyhedra.FullDim(v), eltype(v))`.
"""
function default_library end

_default_type(::Type{T}) where T = T
# See https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/35
_default_type(::Type{AbstractFloat}) = Float64
default_library(::StaticArrays.Size, T::Type) = DefaultLibrary{_default_type(T)}()
default_library(::StaticArrays.Size{(1,)}, T::Type) = IntervalLibrary{_default_type(T)}()
function default_library(d::Int, ::Type{T}) where T
    if d == 1
        return IntervalLibrary{_default_type(T)}()
    else
        return DefaultLibrary{_default_type(T)}()
    end
end

"""
    similar_library(lib::Library, d::FullDim, T::Type)

Returns a library that supports polyhedra of full dimension `T` with coefficient type `T`. If `lib` does not support it, this commonly calls `default_library(d, T)`.
"""
function similar_library end

# Shortcuts
similar_library(p::Polyhedron, d::FullDim, ::Type{T}) where T = similar_library(library(p), d, T)
similar_library(p::Polyhedron, ::Type{T}) where T = similar_library(library(p), T)
similar_library(p::Polyhedron{T}, d::FullDim) where T = similar_library(p, d, T)

"""
    library(p::Polyhedron)

Returns the library used by `p`.
"""
function library end

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
function solver(p::Rep, solver::Union{Nothing, MPB.AbstractMathProgSolver}=default_solver(p))
    if vrepiscomputed(p)
        VRepSolver()
    else
        solver
    end
end

solver(v::VRepresentation, solver::Union{Nothing, MPB.AbstractMathProgSolver}=default_solver(v)) = VRepSolver()
solver(h::HRepresentation, solver::Union{Nothing, MPB.AbstractMathProgSolver}=default_solver(h)) = solver

_promote_reptype(P::Type{<:HRep}, ::Type{<:HRep}) = P
_promote_reptype(P::Type{<:VRep}, ::Type{<:VRep}) = P
# Breaks ambiguity with above two methods
_promote_reptype(P::Type{<:Polyhedron}, ::Type{<:Polyhedron}) = P

function promote_reptype(P1::Type{<:Rep}, P2::Type{<:Rep})
    _promote_reptype(P1, P2)
end
promote_reptype(P1::Type{<:Rep}, P2::Type{<:Rep}, P::Type{<:Rep}...) = promote_reptype(promote_reptype(P1, P2), P...)
promote_reptype(P::Type{<:Rep}) = P

# whether Rep has the constructor Rep(::It...; solver=...)
supportssolver(::Type{<:Rep}) = false

function constructpolyhedron(RepT::Type{<:Rep{T}}, d::FullDim, p::Tuple{Vararg{Rep}}, it::It{T}...) where T
    if supportssolver(RepT)
        solver = default_solver(p...)
        if solver !== nothing
            return RepT(d, it..., solver=solver)
        end
    end
    RepT(d, it...)::RepT # FIXME without this type annotation even convexhull(::PointsHull{2,Int64,Array{Int64,1}}, ::PointsHull{2,Int64,Array{Int64,1}}) is not type stable, why ?
end

function default_similar(p::Tuple{Vararg{Rep}}, d::FullDim, ::Type{T}, it::It{T}...) where T
    # Some types in p may not support `d` or `T` so we call `similar_type` after `promote_reptype`
    RepT = similar_type(promote_reptype(typeof.(p)...), d, T)
    constructpolyhedron(RepT, d, p, it...)
end

"""
    similar(p::Tuple{Vararg{Polyhedra.Rep}}, d::Polyhedra.FullDim, ::Type{T}, it::Polyhedra.It{T}...)

Creates a representation with a type similar to `p` of a polyhedron of full dimension `d`, element type `T` and initialize it with the iterators `it`.
The type of the result will be chosen closer to the type of `p[1]`.
"""
Base.similar(p::Tuple{Vararg{Rep}}, d::FullDim, ::Type{T}, it::It{T}...) where {T} = default_similar(p, d, T, it...)
function promote_coefficient_type(p::Tuple{Vararg{Rep}})
    promote_type(coefficient_type.(p)...)
end
function Base.similar(p::Tuple{Vararg{Rep}}, d::FullDim, it::It...)
    T = promote_coefficient_type(p)
    similar(p, d, T, it...)
end
function Base.similar(p::Tuple{Vararg{Rep}}, it::It...)
    similar(p, FullDim(p[1]), it...)
end
Base.similar(p::Rep, args...) = similar((p,), args...)
