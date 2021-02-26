# # Convex hull and intersection

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Convex hull and intersection.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Convex hull and intersection.ipynb)

# ### Introduction

# These examples illustrate common operations on polyhedra using [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl):
#
# - The convex hull of the union of polytopes
# - The intersection of polytopes
#
# We start by choosing a polyhedral library that will be used for computing the H-representation from the V-representation and vice-versa as well as removing redundant points.
# In these example, we use the default library available in Polyhedra but it can be replaced by any other library listed [here](https://juliapolyhedra.github.io/), e.g. by changing the last two lines below by `import CDDLib` and `lib = CDDLib.Library()` to use [CDDLib](https://github.com/JuliaPolyhedra/CDDLib.jl).

using Test #jl
using Polyhedra
import GLPK
lib = DefaultLibrary{Float64}(GLPK.Optimizer)

# ### Convex hull
#
# The binary convex hull operation between two polyhedra is obtained with the [`convexhull`](https://juliapolyhedra.github.io/Polyhedra.jl/latest/utilities.html#Polyhedra.convexhull) function.
#
# Below we compute the convex hull of the union of two polygons from their V-representation.

P1 = polyhedron(vrep([
    -1.9 -1.7
    -1.8  0.5
     1.7  0.7
     1.9 -0.3
     0.9 -1.1
]), lib)

P2 = polyhedron(vrep([
    -2.5 -1.1
    -0.8  0.8
     0.1  0.9
     1.8 -1.2
     1.3  0.1
]), lib)

Pch = convexhull(P1, P2)

# Note that the convex hull operation is done in the V-representation so no representation conversion is needed for this operation since `P1` and `P2` where constructed from their V-representation:

hrepiscomputed(P1), hrepiscomputed(P2), hrepiscomputed(Pch) #!jl
@test (hrepiscomputed(P1), hrepiscomputed(P2), hrepiscomputed(Pch)) == (false, false, false) #jl

# Let us note that the `convexhull` of a V-representation contains points and rays and represents the convex hull of the points together with the conic hull of the rays. So, `convexhull(P1, P2)` does the union of the vertices:

npoints(Pch) #!jl
@test npoints(Pch) == 10 #jl

# However, if we want to remove the redundant points we can use `removevredundancy!`:

removevredundancy!(Pch)
npoints(Pch) #!jl
@test npoints(Pch) == 8 #jl

# We can plot the polygons and the convex hull of their union using the `plot` function. For further plotting options see the [Plots.jl](http://docs.juliaplots.org/latest/) documentation.
# We can see below the 8 redundant points highlighted with green dots, the two points that are not highlighted are the redundant ones.

using Plots                          #!jl
plot(P1, color="blue", alpha=0.2)    #!jl
plot!(P2, color="red", alpha=0.2)    #!jl
plot!(Pch, color="green", alpha=0.1) #!jl
scatter!(Pch, color="green")         #!jl

# ### Intersection
#
# Intersection of polyhedra is obtained with the [`intersect`](https://juliapolyhedra.github.io/Polyhedra.jl/latest/utilities.html#Base.intersect) function.
#
# Below we compute the intersection of the two polygons from the previous example.

Pint = intersect(P1, P2)
@test nhalfspaces(Pint) == 10 #jl

# While `P1` and `P2` have been constructed from their V-representation, their H-representation has been computed to build the intersection `Pint`.

hrepiscomputed(P1), vrepiscomputed(P1), hrepiscomputed(P2), vrepiscomputed(P2) #!jl
@test (hrepiscomputed(P1), vrepiscomputed(P1), hrepiscomputed(P2), vrepiscomputed(P2)) == (true, true, true, true) #jl

# On the other hand, `Pint` is constructed from its H-representation hence its V-representation has not been computed yet.

hrepiscomputed(Pint), vrepiscomputed(Pint) #!jl
@test (hrepiscomputed(Pint), vrepiscomputed(Pint)) == (true, false) #jl

# We can obtain the number of points in the intersection with `npoints` as follows:

npoints(Pint) #!jl
@test npoints(Pint) == 8 #jl

# Note that this triggers the computation of the V-representation:

hrepiscomputed(Pint), vrepiscomputed(Pint) #!jl
@test (hrepiscomputed(Pint), vrepiscomputed(Pint)) == (true, true) #jl

# We can plot the polygons and their intersection using the `plot` function.

using Plots                            #!jl
plot(P1, color="blue", alpha=0.2)      #!jl
plot!(P2, color="red", alpha=0.2)      #!jl
plot!(Pint, color="yellow", alpha=0.6) #!jl
