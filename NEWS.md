# Polyhedra.jl v0.8 Release Notes

## Load time improvements
- Dependencies on JuMP.jl, RecipesBase.jl, and GeometryBasics.jl were moved to
  weak dependencies on Julia versions supporting package extensions, i.e. v1.9
  and above. On v1.10 this reduces installation time by 15% and load time by
  18% (see [#328]).

## Breaking changes 
The following changes are only breaking on Julia v1.9 and above
- `Polyhedra.Mesh` is now implemented in a package extension requiring
  GeometryBasics.jl. It is sufficient to load your plotting package, i.e.
  Makie.jl or MeshCat.jl, **before** calling `Polyhedra.Mesh`