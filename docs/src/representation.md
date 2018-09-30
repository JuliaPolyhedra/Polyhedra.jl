# Representation

Polyhedra can be described in 2 different ways.

1. H-representation: As the intersection of finitely many halfspaces given by its facets.
2. V-representation: As the convex hull of its vertices + the conic hull of its rays where '+' is the Minkowski sum.

```@docs
HRepresentation
VRepresentation
Representation
```

In `Polyhedra.jl`, those representations are given the respective abstract types `HRepresentation` and `VRepresentation` which are themself subtypes of `Representation`.

These functions can be called on both H-representation and V-representation
```@docs
fulldim
Polyhedra.coefficient_type
```

## H-representation

The fundamental element of an H-representation is the halfspace
```@docs
HalfSpace
```

An H-representation can be created as the intersection of several halfspaces.
For instance, the polytope
```math
\begin{align*}
  x_1 + x_2 &\leq 1 \\
  x_1 - x_2 &\leq 0 \\
  x_1 & \geq 0.
\end{align*}
```
can be created as follows:
```julia
HalfSpace([1, 1], 1) ∩ HalfSpace([1, -1], 0) ∩ HalfSpace([-1, 0], 0)
```

Even if `HalfSpace`s are enough to describe any polyhedron, it is sometimes important to represent the fact that the polyhedron is contained in an affine subspace.
For instance, the simplex
```math
\begin{align*}
  x_1 + x_2 &= 1 \\
  x_1 &\geq 0 \\
  x_2 &\geq 0
\end{align*}
```
is 1-dimensional even if it is defined in a 2-dimensional space.

The fundamental element of an affine subspace is the hyperplane
```@docs
HyperPlane
```

An affine subspace can be created as the intersection of several hyperplanes. For instance
```julia
HyperPlane([1, 1], 1) ∩ HyperPlane([1, 0], 0)
```
represents the 0-dimensional affine subspace only containing the point ``(0, 1)``.

To represent a polyhedron that is not full-dimensional, hyperplanes and halfspaces can be mixed in any order.
For instance, the simplex defined above can be obtained as follows:
```julia
HalfSpace([-1, 0], 0) ∩ HyperPlane([1, 1], 1) ∩ HalfSpace([0, -1], 0)
```

In addition to being created incrementally with intersections, an H-representation can also be created using the `hrep` function
```@docs
hrep
```

### Interface

An H-representation is represented as an intersection halfspaces and hyperplanes. The halfspaces can be obtained with [`halfspaces`](@ref) and the hyperplanes with [`hyperplanes`](@ref).
As an hyperplane ``\langle a, x \rangle = \beta`` is the intersection of the two halfspaces ``\langle a, x \rangle \le \beta`` and ``\langle a, x \rangle \ge \beta``,
even if the H-representation contains hyperplanes, a list of halfspaces whose intersection is the polyhedron can be obtained with [`allhalfspaces`](@ref), which has `nhalfspaces(p) + 2nhyperplanes(p)` elements for an H-representation `p` since each hyperplane is split in two halfspaces.

```@docs
halfspaces
nhalfspaces
hashalfspaces
hyperplanes
nhyperplanes
hashyperplanes
allhalfspaces
nallhalfspaces
hasallhalfspaces
```

## V-representation

The fundamental elements of an V-representation are the points (represented by
`AbstractVector`s and and rays
```@docs
Ray
```

A V-representation can be created as the minkowski sum between a convex hull of points and a conic hull of rays.
For instance, the positive orthant without the simplex defined in the H-representation section can be created as follows:
```julia
convexhull([1, 0], [0, 1]) + conichull([1, 0], [0, 1])
```

The V-representation represents the polyhedron as a minkowski sum of a polytope and a polyhedral cone.
The polytope is represented using a *P-representation* : a convex hull of points.
The polyhedral cone is represented using an *R-representation* : a conic hull of rays.

Even if rays are enough to describe any polyhedral cone, it is sometimes important to represent the fact that the polyhedron contains an affine subspace.
For instance, the polyhedron created with
```julia
convexhull([1, 0], [0, 1]) + conichull([1, 1], [-1, -1])
```
contains the line `[1, 1]`.

The fundamental element of an affine subspace is the line
```@docs
Line
```

An affine subspace can be created as the conic hull/minkownski sum of several lines. For instance
```julia
conichull(Line([1, 0]), Line([0, 1]))
```
represents the full space.

In addition to being created incrementally with convex hull and minkowsky addition, a V-representation can also be created using the `vrep` function
```@docs
vrep
```

### Interface

A P-representation is represented as a convex hull points.
The points can be obtained with [`points`](@ref).

```@docs
points
npoints
haspoints
```

An R-representation is represented as a conic hull of lines and rays.
The rays can be obtained with [`rays`](@ref) and the lines with [`lines`](@ref).
As a line ``r`` is the conic hull of of the two rays ``r`` and ``-r``,
even if the R-representation contains lines, a list of rays whose conic hull is the polyhedral cone can be obtained with [`allrays`](@ref), which has `nrays(R) + 2nlines(R)` elements for an R-representation `R` since each line is split in two rays.

```@docs
rays
nrays
hasrays
lines
nlines
haslines
allrays
nallrays
hasallrays
```
