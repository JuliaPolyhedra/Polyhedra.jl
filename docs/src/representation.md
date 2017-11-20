# Representation

Polyhedra can be described in 2 different ways.

1. H-representation: As the intersection of finitely many halfspaces given by its facets.
2. V-representation: As the convex hull of its vertices + the conic hull of its rays where '+' is the Minkowski sum.

In `Polyhedra.jl`, those representations are given the respective abstract types `HRepresentation` and `VRepresentation` which are themself subtypes of `Representation`.

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

### H-representation interface
```@docs
hreps
nhreps
hashreps
eqs
neqs
haseqs
ineqs
nineqs
hasineqs
hrepiscomputed
```

### V-representation interface
```@docs
vreps
nvreps
hasvreps
rays
nrays
hasrays
points
npoints
haspoints
vrepiscomputed
```
