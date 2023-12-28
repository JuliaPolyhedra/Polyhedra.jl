module Polyhedra

using SparseArrays
using LinearAlgebra
# For `rank` with `Rational{BigInt}` or `BigFloat`,
# it calls `svdvals` for `BigFloat` which is implemented in
# `GenericLinearAlgebra`.
# TODO use exact arithmetic for rank of `Rational{BigInt}`.
#      like RowEchelon.jl (slow) or LU (faster).
import GenericLinearAlgebra

import MutableArithmetics
const MA = MutableArithmetics

import MathOptInterface as MOI
import MathOptInterface.Utilities as MOIU

export Polyhedron

abstract type Library end
abstract type Polyhedron{T} end

export optimizer_with_attributes

coefficient_type(::Union{AbstractVector{T}, Type{<:AbstractVector{T}}}) where T = T
similar_type(::Type{<:Vector}, ::Int, ::Type{T}) where T = Vector{T}
similar_type(::Type{SparseVector{S, IT}}, ::Int, ::Type{T}) where {S, IT, T} = SparseVector{T, IT}

emptymatrix(::Type{MT}, m, n) where {MT<:AbstractMatrix} = MT(undef, m, n)
emptymatrix(::Type{SparseMatrixCSC{T, Int}}, m, n) where T = spzeros(T, m, n)
similar_type(::Type{<:Matrix}, ::Type{T}) where T = Matrix{T}
similar_type(::Type{SparseMatrixCSC{S, I}}, ::Type{T}) where {S, I, T} = SparseMatrixCSC{T, I}

# Interface/Definitions
include("dimension.jl")
include("elements.jl")
include("comp.jl")
include("representation.jl")
include("indices.jl")
include("incidence.jl")
include("triangulation.jl")
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
vectortype(::Type{StaticArrays.SArray{Tuple{N, M}, T, 2, NM}}) where {N, M, T, NM} = StaticArrays.SVector{M, T}

hmatrixtype(RepT::Type{<:HRep}, T::Type) = matrixtype(similar_type(hvectortype(RepT), T))
vmatrixtype(RepT::Type{<:VRep}, T::Type) = matrixtype(similar_type(vvectortype(RepT), T))
matrixtype(::Type{<:AbstractVector{T}}) where T = Matrix{T}
matrixtype(::Type{<:AbstractSparseVector{T}}) where T = SparseMatrixCSC{T, Int}

function similar_type(::Type{ET}, ::Type{Tout}) where {Tout, ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, FullDim(ET), Tout)
end
function similar_type(::Type{ET}, d::FullDim) where {ET<:Union{HRepElement, VRepElement, Rep}}
    similar_type(ET, d, coefficient_type(ET))
end

# Operations
include("repop.jl")
include("center.jl")
include("repelemop.jl")
include("aff.jl")
include("linearity.jl")
include("redundancy.jl")
include("projection.jl")
include("extended.jl")

# Implementations of representations
include("vecrep.jl")
include("mixedrep.jl")
include("lphrep.jl")
include("matrep.jl")
include("liftedrep.jl")
include("doubledescription.jl") # FIXME move it after projection.jl once it stops depending on LiftedRep
include("interval.jl") # 1D polyhedron
include("planar.jl")
include("defaultlibrary.jl")

include("fulldim.jl")

# Optimization
include("opt.jl")
include("polyhedra_to_lp_bridge.jl")
include("vrep_optimizer.jl")
include("default.jl")
include("projection_opt.jl")

# Visualization
include("show.jl")

"""
    Mesh(p::Polyhedron)

Returns wrapper of a polyhedron suitable for plotting with MeshCat.jl and Makie.jl.

!!! note "Extension in Julia 1.9 and above"
    Although we require `using GeometryBasics` to use this function in Julia 1.9 and above,
    in most use cases this extension dependency is loaded by the plotting package and no
    further action is required.
"""
function Mesh end

if !isdefined(Base, :get_extension)
    include("../ext/PolyhedraJuMPExt.jl")
    include("../ext/PolyhedraRecipesBaseExt.jl")
    include("../ext/PolyhedraGeometryBasicsExt.jl")
end

end # module
