__precompile__()

module Polyhedra

using FixedSizeArrays, GeometryTypes

export PolyhedraLibrary, Polyhedron, getlibrary, getlibraryfor

abstract PolyhedraLibrary
abstract Polyhedron{N,T} <: GeometryPrimitive{N,T}

import Base.intersect, Base.+, Base.*, Base.isempty, Base.copy, Base.push!, Base.length, Base.eltype, Base.start, Base.done, Base.next

# Definitions
include("elements.jl")
include("mycomp.jl")
include("representation.jl")
include("repop.jl")
include("repelemop.jl")
include("operations.jl")
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
include("decompose.jl")

end # module
