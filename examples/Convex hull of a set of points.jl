# # Convex hull of a set of points

# ## Convex hull in a plane

# This examples shows how to find the convex hull in the context of a plane.
# First we have to create an object representing the vertices. We have multiple ways of doing this

using Polyhedra

v = convexhull([0, 0], [0, 1], [1, 0], [0.1, 0.1]) # list of points
v = vrep([[0, 0], [0, 1], [1, 0], [0.1, 0.1]]) # vector of points
# number of points × dimension matrix
x = [0,0,1,0.1]
y = [0,1,0,0.1]
v = vrep([x y])

# Then we can compute the hull of these points using the planar_hull function

Polyhedra.planar_hull(v)

# ## Convex hull in higher dimension

# In higher dimension, we can do it with a linear programming solver implementing the [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl), e.g.,

import GLPK

removevredundancy(v, GLPK.Optimizer)

# We can also use any Polyhedral library implementing the interface of this package. If we don't specify any library, it falls back to a default one implementing on this package which will use the `planar_hull` if the dimension is 2 (so it's equivalent to the first approach shown above):

p = polyhedron(v)
removevredundancy!(p)
p

# We can also specify a library:

using CDDLib

p = polyhedron(v, CDDLib.Library())
removevredundancy!(p)
p
