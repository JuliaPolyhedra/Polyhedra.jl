# Polyhedra --- Manipulation of Polyhedra in Julia

[Polyhedra](https://github.com/JuliaPolyhedra/Polyhedra.jl) is a package for polyhedra manipulations in Julia.
It provides an unified interface for Polyhedra Manipulation Libraries such as [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl), [LRSLib.jl](https://github.com/JuliaPolyhedra/LRSLib.jl) and [QHull](https://github.com/davidavdav/QHull.jl).

Polyhedra can either be represented by a set of linear inequalities or by vertices and rays.
In the first case, the points of the polyhedron are the points that satisfies all the inequalities
and in the second case they are the points that can be expressed as a convex combination of the vertices plus a conic combination of the rays.
The manipulations that `Polyhedra` can perform include

* Projection: Projection of a polyhedron on a lower dimensional space, e.g. Fourier-Motzkin elimination.
* Changing the Representation

  * Vertex enumeration problem: Computing the extremal vertices and rays from an inequality representation
  * Convex hull problem: Computing a set of linear inequalities describing the polyhedron from a vertex/ray representation

* Removal of redundant inequalities or redundant vertices/rays.
* Plotting of 2D polyhedra using [Plots](https://github.com/JuliaPlots/Plots.jl)
* Decomposition of 3D polyhedra into vertices and triangular faces,
  enabling easy visualization of 3D polyhedra using
  [DrakeVisualizer](https://github.com/rdeits/DrakeVisualizer.jl) or
  [GLVisualize](https://github.com/JuliaGL/GLVisualize.jl).

Depending on the library, these manipulation can either be in floating point or exact rational arithmetic.

Each operation has a default fallback implementation which is used in case the library does not support it.
Polyhedra also includes a default library which does not implement anything, hence using every fallback.

Polyhedra remains under active development, and we welcome your feedback, suggestions, and bug reports.

## Installing Polyhedra

If you are familiar with Julia you can get started quickly by using the
package manager to install Polyhedra
```julia
julia> Pkg.add("Polyhedra")
```

And a Polyhedra Manipulation Library, e.g.
```julia
julia> Pkg.add("CDDLib")
```

## Contents
```@contents
Pages = ["installation.md", "representation.md", "polyhedron.md", "plot.md", "redundancy.md", "projection.md", "optimization.md", "utilities.md"]
Depth = 2
```
