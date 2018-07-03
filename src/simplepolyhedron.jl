export SimplePolyhedraLibrary, SimplePolyhedron

"""
    SimplePolyhedraLibrary{T}

Default library for polyhedra of dimension larger than 1 ([`IntervalLibrary`](@ref) is the default for polyhedra of dimension 1).
The library implements the bare minimum and uses the fallback implementation for all operations.
"""
struct SimplePolyhedraLibrary{T} <: PolyhedraLibrary
    solver::Union{Nothing, MPB.AbstractMathProgSolver}
    function SimplePolyhedraLibrary{T}(solver=nothing) where T
        new{T}(solver)
    end
end

similar_library(lib::SimplePolyhedraLibrary, d::FullDim, ::Type{T}) where T = default_library(d, T) # default_library allows to fallback to Interval if d is FullDim{1}

mutable struct SimplePolyhedron{N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}} <: Polyhedron{N, T}
    hrep::Nullable{HRepT}
    vrep::Nullable{VRepT}
    solver::Union{Nothing, MPB.AbstractMathProgSolver}
    function SimplePolyhedron{N, T, HRepT, VRepT}(hrep::Union{HRepT, Void}, vrep::Union{VRepT, Void}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
        new{N, T, HRepT, VRepT}(hrep, vrep, solver)
    end
end
function SimplePolyhedron{N, T, HRepT, VRepT}(hrep::HRepT, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
    SimplePolyhedron{N, T, HRepT, VRepT}(hrep, nothing, solver)
end
function SimplePolyhedron{N, T, HRepT, VRepT}(vrep::VRepT, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
    SimplePolyhedron{N, T, HRepT, VRepT}(nothing, vrep, solver)
end
function SimplePolyhedron{N, T, HRepT, VRepT}(hrep::HRepresentation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
    SimplePolyhedron{N, T, HRepT, VRepT}(HRepT(hrep), solver)
end
function SimplePolyhedron{N, T, HRepT, VRepT}(vrep::VRepresentation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T, HRepT<:HRepresentation{N, T}, VRepT<:VRepresentation{N, T}}
    SimplePolyhedron{N, T, HRepT, VRepT}(VRepT(vrep), solver)
end

library(::Union{SimplePolyhedron{N, T}, Type{<:SimplePolyhedron{N, T}}}) where {N, T} = SimplePolyhedraLibrary{T}()
default_solver(p::SimplePolyhedron) = p.solver
supportssolver(::Type{<:SimplePolyhedron}) = true

hvectortype(::Type{<:SimplePolyhedron{N, T, HRepT}}) where {N, T, HRepT} = hvectortype(HRepT)
vvectortype(::Type{SimplePolyhedron{N, T, HRepT, VRepT}}) where {N, T, HRepT, VRepT} = vvectortype(VRepT)

similar_type(::Type{<:SimplePolyhedron{M, S, HRepT, VRepT}}, d::FullDim{N}, ::Type{T}) where {M, S, HRepT, VRepT, N, T} = SimplePolyhedron{N, T, similar_type(HRepT, d, T), similar_type(VRepT, d, T)}

function SimplePolyhedron{N, T, HRepT, VRepT}(hits::HIt...; solver=nothing) where {N, T, HRepT, VRepT}
    SimplePolyhedron{N, T, HRepT, VRepT}(HRepT(hits...), solver)
end
function SimplePolyhedron{N, T, HRepT, VRepT}(vits::VIt...; solver=nothing) where {N, T, HRepT, VRepT}
    SimplePolyhedron{N, T, HRepT, VRepT}(VRepT(vits...), solver)
end

# Need fulltype in case the use does `intersect!` with another element
SimplePolyhedron{N, T}(rep::Representation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T} = SimplePolyhedron{N, T}(MultivariatePolynomials.changecoefficienttype(rep, T), solver)
function SimplePolyhedron{N, T}(rep::HRepresentation{N, T}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T}
    HRepT = fulltype(typeof(rep))
    VRepT = dualtype(HRepT)
    SimplePolyhedron{N, T, HRepT, VRepT}(rep, solver)
end
function SimplePolyhedron{N, T}(rep::VRepresentation{N, T}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {N, T}
    VRepT = fulltype(typeof(rep))
    HRepT = dualtype(VRepT)
    SimplePolyhedron{N, T, HRepT, VRepT}(rep, solver)
end

function polyhedron(rep::Representation{N}, lib::SimplePolyhedraLibrary{T}) where {N, T}
    SimplePolyhedron{N, polytypefor(T)}(rep, lib.solver)
end

function Base.copy(p::SimplePolyhedron{N, T}) where {N, T}
    if !isnull(p.hrep)
        SimplePolyhedron{N, T}(get(p.hrep), p.solver)
    else
        SimplePolyhedron{N, T}(get(p.vrep), p.solver)
    end
end

hrepiscomputed(p::SimplePolyhedron) = !isnull(p.hrep)
function computehrep!(p::SimplePolyhedron)
    # vrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.hrep = doubledescription(get(p.vrep))
end
function hrep(p::SimplePolyhedron)
    if !hrepiscomputed(p)
        computehrep!(p)
    end
    get(p.hrep)
end
vrepiscomputed(p::SimplePolyhedron) = !isnull(p.vrep)
function computevrep!(p::SimplePolyhedron)
    # hrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.vrep = doubledescription(get(p.hrep))
end
function vrep(p::SimplePolyhedron)
    if !vrepiscomputed(p)
        computevrep!(p)
    end
    get(p.vrep)
end

function sethrep!(p::SimplePolyhedron, h::HRepresentation)
    p.hrep = h
end
function setvrep!(p::SimplePolyhedron, v::VRepresentation)
    p.vrep = v
end
function resethrep!(p::SimplePolyhedron, h::HRepresentation)
    p.hrep = h
    p.vrep = nothing
end
function resetvrep!(p::SimplePolyhedron, v::VRepresentation)
    p.vrep = v
    p.hrep = nothing
end
