# Representation

Polyhedra can be described in 2 different ways.

1. H-representation: As the intersection of finitely many halfspaces given by its facets.
2. V-representation: As the convex hull of its vertices + the conic hull of its rays where '+' is the Minkowski sum.

In `Polyhedra.jl`, those representations are given the respective abstract types `HRepresentation` and `VRepresentation` which are themself subtypes of `Representation`.

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
For instance, the simplex:
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

An H-representation is represented as an intersection halfspaces and hyperplanes. The halfspaces can be obtained with `halfspaces` and the hyperplanes with `hyperplane`.
As an hyperplane ``\langle a, x \rangle = \beta`` is the intersection of the two halfspaces ``\\\\langle a, x \\\\rangle \\\\le \\\\beta`` and ``\\langle a, x \\rangle \\ge \\beta``,
even if the H-representation contains hyperplanes, a list of halfspaces whose intersection is the polyhedron can be obtained with `allhalfspaces`, which has `nhalfspaces(p) + 2nhyperplanes(p)` elements for an H-representation `p` since each hyperplane is split in two halfspaces.

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
hrepiscomputed
```

## V-representation

For instance, consider the 2-dimensional polyhedron described by the following H-representation:
```math
\begin{align*}
  x_1 + x_2 &\leq 1 \\
  x_1 - x_2 &\leq 0 \\
  x_1 & \geq 0.
\end{align*}
```

This set of inequalities can be written in the matrix form ``Ax \leq b`` where
```math
A = \begin{pmatrix}1 & 1\\1 & -1\\-1 & 0\end{pmatrix}, b = \begin{pmatrix}1\\0\\0\end{pmatrix}.
```

Let's create this H-representation using the concrete subtype `SimpleHRepresentation` of the abstract type `HRepresentation`.
```julia
julia> using Polyhedra
julia> A = [1 1;1 -1;-1 0]
julia> b = [1,0,0]
julia> hrep = SimpleHRepresentation(A, b)
julia> typeof(hrep)
Polyhedra.SimpleHRepresentation{2,Int64}
```

This polyhedron has three vertices: ``(0,0)``, ``(0,1)`` and ``(0.5,0.5)``.
We can create this V-representation using the concrete subtype `SimpleVRepresentation` of the abstract type `VRepresentation`.
Because ``0.5`` is fractional, have two choices: either use exact rational arithemtic
```julia
julia> V = [0 0; 0 1; 1//2 1//2]
julia> vrep = SimpleVRepresentation(V)
julia> typeof(vrep)
Polyhedra.SimpleVRepresentation{2,Rational{Int64}}
```

or use floating point arithmetic
```julia
julia> Vf = [0 0; 0 1; 1/2 1/2]
julia> vrepf = SimpleVRepresentation(Vf)
julia> typeof(vrepf)
Polyhedra.SimpleVRepresentation{2,Float64}
```

## Representation interface
These functions can be called on both H-representation and V-representation
```@docs
fulldim
```

### V-representation interface
```@docs
points
npoints
haspoints
sympoints
nsympoints
hassympoints
allpoints
nallpoints
hasallpoints
rays
nrays
hasrays
lines
nlines
haslines
allrays
nallrays
hasallrays
vrepiscomputed
```
