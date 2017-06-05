__precompile__()

module Polyhedra

using GeometryTypes
using StaticArrays.FixedSizeArrays: FixedVector

export PolyhedraLibrary, Polyhedron, getlibrary, getlibraryfor

abstract type PolyhedraLibrary end
abstract type Polyhedron{N,T} <: GeometryPrimitive{N,T} end

import Base.intersect, Base.==, Base.+, Base.*, Base.isempty, Base.copy, Base.push!, Base.length, Base.eltype, Base.start, Base.done, Base.next

# Definitions
include("elements.jl")
include("mycomp.jl")
include("representation.jl")
include("repop.jl")
include("repelemop.jl")
include("operations.jl")
include("aff.jl")
include("redundancy.jl")
include("projection.jl")

# Implementations
include("lphrep.jl")
include("jump.jl")
include("simplerep.jl")
include("liftedrep.jl")
include("doubledescription.jl")
include("simplepolyhedron.jl")

# Optimization
importall MathProgBase.SolverInterface
include("opt.jl")
include("lpqp_to_polyhedra.jl")
include("polyhedra_to_lpqp.jl")
include("simplevrepsolver.jl")
include("default.jl")

# Visualization
include("show.jl")
include("recipe.jl")
include("decompose.jl")

end # module
