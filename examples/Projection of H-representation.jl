# # Projection of H-representation

# This examples shows how to compute the projection of a H-representation.
# We use the [project1.ine example from cddlib](https://github.com/cddlib/cddlib/blob/master/examples/project1.ine).

using Polyhedra

# ## Projection with default library

h = LiftedHRepresentation{Float64}([
    1  0  0  0  1  0  0
    1  0  0  0  0  1  0
    1  0  0  0  0  0  1
    1  0  0  0 -1  0  0
    1  0  0  0  0 -1  0
    1  0  0  0  0  0 -1
    1  1  0  0 -1  0  0
    1  0  1  0  0 -1  0
    1  0  0  1  0  0 -1
    1 -1  0  0  1  0  0
    1  0 -1  0  0  1  0
    1  0  0 -1  0  0  1
    2  1  1  1 -1 -1 -1
    2 -1  1  1  1 -1 -1
    2  1 -1  1 -1  1 -1
    2  1  1 -1 -1 -1  1
    2 -1 -1  1  1  1 -1
    2  1 -1 -1 -1  1  1
    2 -1  1 -1  1 -1  1
    2 -1 -1 -1  1  1  1
])
p = polyhedron(convert(Polyhedra.Intersection{Float64,Vector{Float64},Int}, h))

# The projection is going to first compute the V-representation,
# and then project this V-representation. Only the V-representation
# of the projected polyhedron will be known.

project(p, 1:3)

# ## CDDLib
#
# With CDDLib, you the default projection when the eliminated variables is
# not only the last dimension is [`BlockElimination`](@ref):
# For this method, the H-representation of the projected set is obtained:

import CDDLib
p = polyhedron(h, CDDLib.Library())

# [`FourierMotzkin`](@ref) can be used instead as follows.
# This method also obtains the H-representation of the projected set.

project(p, 1:3, FourierMotzkin())

# To compute the projection by first computing the V-representation
# and then projecting it like it was done for the default library,
# use [`ProjectGenerators`](@ref).

project(p, 1:3, ProjectGenerators())
