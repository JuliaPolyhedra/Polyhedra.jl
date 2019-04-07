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
function default_type(d::Integer, ::Type{T}) where T
    if isone(d)
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
function default_library(d::Integer, ::Type{T}) where T
    if isone(d)
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
default_solver(p::Rep; kws...) = nothing
function default_solver(p::Rep, ps::Rep...; kws...)
    s = default_solver(p; kws...)
    if s === nothing
        return default_solver(ps...; kws...)
    else
        return s
    end
end

"""
    linear_objective_solver(p::Rep, solver::Union{Nothing, JuMP.OptimizerFactory}=default_solver(p))

Return the solver to use for optimizing a linear objective over the polyhedron `p`, i.e.
```julia
model = Model(solver)
x = @variable(model, [1:fulldim(p)])
@constraint(model, x in p)
@objective(model, c ⋅ x)
```
for some vector `c`.

By default, if the V-representation of `p` has been computed, it returns
`VRepOptimizer()`, otherwise, it returns `solver`.

If the problem has constraints different to `x in p`, use `default_solver(p)` instead
as the fact that the V-representation of `p` has been computed does not help.
"""
function linear_objective_solver(p::Rep{T}, solver::SolverOrNot=default_solver(p)) where T
    if vrepiscomputed(p)
        return with_optimizer(VRepOptimizer{T})
    else
        return solver
    end
end

linear_objective_solver(v::VRepresentation{T}, solver::SolverOrNot=default_solver(v)) where {T} = with_optimizer(VRepOptimizer{T})
linear_objective_solver(h::HRepresentation, solver::SolverOrNot=default_solver(h)) = solver

_promote_reptype(P::Type{<:HRep}, ::Type{<:HRep}) = P
_promote_reptype(P::Type{<:VRep}, ::Type{<:VRep}) = P
# Breaks ambiguity with above two methods
_promote_reptype(P::Type{<:Polyhedron}, ::Type{<:Polyhedron}) = P

function promote_reptype(P1::Type{<:Rep}, P2::Type{<:Rep})
    return _promote_reptype(P1, P2)
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
    return RepT(d, it...)::RepT # FIXME without this type annotation even convexhull(::PointsHull{2,Int64,Array{Int64,1}}, ::PointsHull{2,Int64,Array{Int64,1}}) is not type stable, why ?
end

function default_similar(p::Tuple{Vararg{Rep}}, d::FullDim, ::Type{T}, it::It{T}...) where T
    # Some types in p may not support `d` or `T` so we call `similar_type` after `promote_reptype`
    RepT = similar_type(promote_reptype(typeof.(p)...), d, T)
    return constructpolyhedron(RepT, d, p, it...)
end

"""
    similar(p::Tuple{Vararg{Polyhedra.Rep}}, d::Polyhedra.FullDim, ::Type{T}, it::Polyhedra.It{T}...)

Creates a representation with a type similar to `p` of a polyhedron of full dimension `d`, element type `T` and initialize it with the iterators `it`.
The type of the result will be chosen closer to the type of `p[1]`.
"""
Base.similar(p::Tuple{Vararg{Rep}}, d::FullDim, ::Type{T}, it::It{T}...) where {T} = default_similar(p, d, T, it...)
function promote_coefficient_type(p::Tuple{Vararg{Rep}})
    return promote_type(coefficient_type.(p)...)
end
function Base.similar(p::Tuple{Vararg{Rep}}, d::FullDim, it::It...)
    T = promote_coefficient_type(p)
    return similar(p, d, T, it...)
end
function Base.similar(p::Tuple{Vararg{Rep}}, it::It...)
    return similar(p, FullDim(p[1]), it...)
end
Base.similar(p::Rep, args...) = similar((p,), args...)
