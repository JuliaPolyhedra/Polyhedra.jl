__precompile__()

module Polyhedra

using Compat
using Compat.SparseArrays
using Compat.LinearAlgebra

using MultivariatePolynomials

export PolyhedraLibrary, Polyhedron

abstract type PolyhedraLibrary end
abstract type Polyhedron{T} end

import MathProgBase
const MPB = MathProgBase
const MPBSI = MPB.SolverInterface

include("dimension.jl")

MultivariatePolynomials.coefficienttype(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::Int, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Int, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}

emptymatrix(::Type{MT}, m, n) where {MT<:AbstractMatrix} = MT(undef, m, n)
emptymatrix(::Type{SparseMatrixCSC{T, Int}}, m, n) where T = spzeros(T, m, n)
similar_type(::Type{<:Matrix}, ::Type{T}) where T = Matrix{T}
similar_type(::Type{SparseMatrixCSC{S, I}}, ::Type{T}) where {S, I, T} = SparseMatrixCSC{T, I}

# Interface/Definitions
include("elements.jl")
include("comp.jl")
include("representation.jl")
include("indices.jl")
include("incidence.jl")
include("iterators.jl")
include("polyhedron.jl")

vvectortype(rep::Rep) = vvectortype(typeof(rep))
hvectortype(rep::Rep) = hvectortype(typeof(rep))
# TODO Only define vectortype for the type for Polyhedron v0.4
function vectortype(RepT::Type{<:Polyhedron})
    @assert hvectortype(RepT) == vvectortype(RepT)
    hvectortype(RepT)
end
vectortype(HRepT::Type{<:HRepresentation}) = hvectortype(HRepT)
vectortype(VRepT::Type{<:VRepresentation}) = vvectortype(VRepT)
vectortype(p::Rep) = vectortype(typeof(p))

vectortype(::Type{<:AbstractSparseArray{T}}) where T = SparseVector{T, Int}
vectortype(::Type{<:AbstractMatrix{T}}) where T = Vector{T}

hmatrixtype(RepT::Type{<:HRep}, T::Type) = matrixtype(similar_type(hvectortype(RepT), T))
vmatrixtype(RepT::Type{<:VRep}, T::Type) = matrixtype(similar_type(vvectortype(RepT), T))
matrixtype(::Type{<:AbstractVector{T}}) where T = Matrix{T}
matrixtype(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, d, coefficienttype(ET))
end

# Operations
include("repop.jl")
include("center.jl")
include("repelemop.jl")
include("aff.jl")
include("redundancy.jl")
include("projection.jl")

# Implementations of representations
include("vecrep.jl")
include("mixedrep.jl")
include("lphrep.jl")
include("jump.jl")
include("matrep.jl")
include("liftedrep.jl")
include("doubledescription.jl") # FIXME move it after projection.jl once it stops depending on LiftedRep
include("interval.jl") # 1D polyhedron
include("simplepolyhedron.jl")

# -1 is the dimension of an empty polyhedron, here it is used as the
# *full* dimension of a polyhedron with no element
fulldim_rec() = -1
function fulldim_rec(rep::Rep{T}, its::Union{Rep{T}, It{T}}...) where T
    N = fulldim(rep)
    if N == -1
        return fulldim_rec(its...)
    else
        return N
    end
end
function fulldim_rec(it::It{T}, its::Union{Rep{T}, It{T}}...) where T
    if isempty(it)
        return fulldim_rec(its...)
    else
        return fulldim(first(it))
    end
end
function FullDim(::Union{HyperPlanesIntersection{T, AT},
                         LinesHull{T, AT},
                         VEmptySpace{T, AT},
                         Intersection{T, AT},
                         PointsHull{T, AT},
                         RaysHull{T, AT},
                         Hull{T, AT}}) where {T, AT <: StaticArrays.SVector}
    return FullDim(AT)
end



# Optimization
include("opt.jl")
include("lpqp_to_polyhedra.jl")
include("polyhedra_to_lpqp.jl")
include("vrepsolver.jl")
include("default.jl")

# Visualization
include("show.jl")
include("recipe.jl")
include("decompose.jl")

end # module
