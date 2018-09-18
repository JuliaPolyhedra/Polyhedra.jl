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

mutable struct SimplePolyhedron{T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}} <: Polyhedron{T}
    hrep::Union{HRepT, Nothing}
    vrep::Union{VRepT, Nothing}
    solver::Union{Nothing, MPB.AbstractMathProgSolver}
    function SimplePolyhedron{T, HRepT, VRepT}(hrep::Union{HRepT, Nothing}, vrep::Union{VRepT, Nothing}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}}
        new{T, HRepT, VRepT}(hrep, vrep, solver)
    end
end
function SimplePolyhedron{T, HRepT, VRepT}(hrep::HRepT, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}}
    SimplePolyhedron{T, HRepT, VRepT}(hrep, nothing, solver)
end
function SimplePolyhedron{T, HRepT, VRepT}(vrep::VRepT, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}}
    SimplePolyhedron{T, HRepT, VRepT}(nothing, vrep, solver)
end
function SimplePolyhedron{T, HRepT, VRepT}(hrep::HRepresentation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}}
    SimplePolyhedron{T, HRepT, VRepT}(HRepT(hrep), solver)
end
function SimplePolyhedron{T, HRepT, VRepT}(vrep::VRepresentation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T, HRepT<:HRepresentation{T}, VRepT<:VRepresentation{T}}
    SimplePolyhedron{T, HRepT, VRepT}(VRepT(vrep), solver)
end

library(::Union{SimplePolyhedron{T}, Type{<:SimplePolyhedron{T}}}) where {T} = SimplePolyhedraLibrary{T}()
default_solver(p::SimplePolyhedron) = p.solver
supportssolver(::Type{<:SimplePolyhedron}) = true

hvectortype(::Type{<:SimplePolyhedron{T, HRepT}}) where {T, HRepT} = hvectortype(HRepT)
vvectortype(::Type{SimplePolyhedron{T, HRepT, VRepT}}) where {T, HRepT, VRepT} = vvectortype(VRepT)

similar_type(::Type{<:SimplePolyhedron{M, S, HRepT, VRepT}}, d::FullDim{N}, ::Type{T}) where {M, S, HRepT, VRepT, T} = SimplePolyhedron{T, similar_type(HRepT, d, T), similar_type(VRepT, d, T)}

function SimplePolyhedron{T, HRepT, VRepT}(hits::HIt...; solver=nothing) where {T, HRepT, VRepT}
    SimplePolyhedron{T, HRepT, VRepT}(HRepT(hits...), solver)
end
function SimplePolyhedron{T, HRepT, VRepT}(vits::VIt...; solver=nothing) where {T, HRepT, VRepT}
    SimplePolyhedron{T, HRepT, VRepT}(VRepT(vits...), solver)
end

# Need fulltype in case the use does `intersect!` with another element
SimplePolyhedron{T}(rep::Representation{N}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T} = SimplePolyhedron{T}(MultivariatePolynomials.changecoefficienttype(rep, T), solver)
function SimplePolyhedron{T}(rep::HRepresentation{T}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T}
    HRepT = fulltype(typeof(rep))
    VRepT = dualtype(HRepT)
    SimplePolyhedron{T, HRepT, VRepT}(rep, solver)
end
function SimplePolyhedron{T}(rep::VRepresentation{T}, solver::Union{Nothing, MPB.AbstractMathProgSolver}) where {T}
    VRepT = fulltype(typeof(rep))
    HRepT = dualtype(VRepT)
    SimplePolyhedron{T, HRepT, VRepT}(rep, solver)
end

function polyhedron(rep::Representation{N}, lib::SimplePolyhedraLibrary{T}) where {T}
    SimplePolyhedron{polytypefor(T)}(rep, lib.solver)
end

function Base.copy(p::SimplePolyhedron{T}) where {T}
    if p.hrep !== nothing
        SimplePolyhedron{N, T}(p.hrep, p.solver)
    else
        SimplePolyhedron{N, T}(p.vrep, p.solver)
    end
end

hrepiscomputed(p::SimplePolyhedron) = p.hrep !== nothing
function computehrep!(p::SimplePolyhedron)
    # vrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.hrep = doubledescription(p.vrep)
end
function hrep(p::SimplePolyhedron)
    if !hrepiscomputed(p)
        computehrep!(p)
    end
    return p.hrep
end
vrepiscomputed(p::SimplePolyhedron) = p.vrep !== nothing
function computevrep!(p::SimplePolyhedron)
    # hrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.vrep = doubledescription(p.hrep)
end
function vrep(p::SimplePolyhedron)
    if !vrepiscomputed(p)
        computevrep!(p)
    end
    return p.vrep
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
