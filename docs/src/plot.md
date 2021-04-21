# Plot

Polyhedra contains utilities to visualize either a 2-dimensional or a 3-dimensional
polyhedron, see [Polyhedron](@ref) for how to construct a polyhedron, e.g. from its H- or V-representation.

## 2D plotting with Plots

A 2-dimensional polytope, i.e. *bounded* polyhedron, can be visualized with [Plots](https://github.com/JuliaPlots/Plots.jl).
Suppose for instance that we want to visualize the polyhedron having the following H-representation:
```@example plots2
using Polyhedra
h = HalfSpace([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0)
```

The H-representation cannot be given to Plots directly, it first need to be transformed into a polyhedron:
```@example plots2
p = polyhedron(h)
```

The polyhedron can be given to Plots as follows.
We use `ratio=:equal` so that the horizontal and vertical axis have the same scale.
```@example plots2
using Plots
plot(p, ratio=:equal)
```

See [Polyhedral Function](https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/master/examples/Polyhedral%20Function.ipynb) and [3D Plotting a projection of the 4D permutahedron](https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/master/examples/3D%20Plotting%20a%20projection%20of%20the%204D%20permutahedron.ipynb) for example notebooks.

## 3D plotting with Plots

A 3-dimensional polyhedron can be visualized with either [MeshCat](https://github.com/rdeits/MeshCat.jl) or [Makie](https://github.com/JuliaPlots/Makie.jl).
Unbounded polyhedron are supported by truncating the polyhedron into a polytope and not triangularizing the faces in the directions of unbounded rays.

Suppose for instance that we want to visualize the polyhedron having the following H-representation:
```jldoctest plots3
julia> using Polyhedra

julia> v = convexhull([0, 0, 0]) + conichull([1, 0, 0], [0, 1, 0], [0, 0, 1])
V-representation Polyhedra.Hull{Int64, Vector{Int64}, Int64}:
1-element iterator of Vector{Int64}:
 [0, 0, 0],
3-element iterator of Ray{Int64, Vector{Int64}}:
 Ray([1, 0, 0])
 Ray([0, 1, 0])
 Ray([0, 0, 1])
```

The V-representation cannot be given to [MeshCat](https://github.com/rdeits/MeshCat.jl) or [Makie](https://github.com/JuliaPlots/Makie.jl) directly, it first need to be transformed into a polyhedron:
```jldoctest plots3
julia> p = polyhedron(v)
Polyhedron DefaultPolyhedron{Rational{BigInt}, Polyhedra.Intersection{Rational{BigInt}, Vector{Rational{BigInt}}, Int64}, Polyhedra.Hull{Rational{BigInt}, Vector{Rational{BigInt}}, Int64}}:
1-element iterator of Vector{Rational{BigInt}}:
 Rational{BigInt}[0//1, 0//1, 0//1],
3-element iterator of Ray{Rational{BigInt}, Vector{Rational{BigInt}}}:
 Ray(Rational{BigInt}[1//1, 0//1, 0//1])
 Ray(Rational{BigInt}[0//1, 1//1, 0//1])
 Ray(Rational{BigInt}[0//1, 0//1, 1//1])
```

Then, we need to create a mess from the polyhedron:
```jldoctest plots3
julia> m = Polyhedra.Mesh(p)
Polyhedra.Mesh{3, Rational{BigInt}, DefaultPolyhedron{Rational{BigInt}, Polyhedra.Intersection{Rational{BigInt}, Vector{Rational{BigInt}}, Int64}, Polyhedra.Hull{Rational{BigInt}, Vector{Rational{BigInt}}, Int64}}}(convexhull([0//1, 0//1, 0//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1, 0//1]), Ray(Rational{BigInt}[0//1, 1//1, 0//1]), Ray(Rational{BigInt}[0//1, 0//1, 1//1])), nothing, nothing, nothing)
```

```@docs
Polyhedra.Mesh
```

The polyhedron can be plotted with [MeshCat](https://github.com/rdeits/MeshCat.jl) as follows
```julia
julia> using MeshCat

julia> vis = Visualizer()

julia> setobject!(vis, m)

julia> open(vis)
```
To plot it in a notebook, replace `open(vis)` with `IJuliaCell(vis)`.

To plot it with [Makie](https://github.com/JuliaPlots/Makie.jl) instead, you can use for instance `mesh` or `wireframe`.
```julia
julia> import Makie

julia> Makie.mesh(m, color=:blue)

julia> Makie.wireframe(m)
```

See [3D Plotting a projection of the 4D permutahedron](https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/master/examples/3D%20Plotting%20a%20projection%20of%20the%204D%20permutahedron.ipynb) for an example notebook.
