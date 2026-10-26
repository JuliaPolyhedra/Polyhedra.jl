# Polyhedra.jl upcoming release notes

## Breaking changes
- The dimension query `Polyhedra.FullDim(p)` is now
  `Polyhedra.typed_fulldim(p)`. Implementations of the Polyhedra interface
  should extend `typed_fulldim`; `FullDim` remains the dimension type alias.
- The undefined `start*`, `done*`, and `next*` iterator names are no longer
  exported. Use Julia's `iterate` interface instead.

# Polyhedra.jl v0.8 Release Notes

## Load time improvements
- Dependencies on JuMP.jl, RecipesBase.jl, and GeometryBasics.jl were moved to
  weak dependencies on Julia versions supporting package extensions, i.e. v1.9
  and above. On v1.10 this reduces installation time by 15% and load time by
  18% (see [#328]).

## Breaking changes
- `JuMP.optimizer_with_attributes` is no longer exported. Call it from JuMP.jl instead.
- The following change is only breaking on Julia v1.9 and above:
  `Polyhedra.Mesh` is now implemented in a package extension requiring
  GeometryBasics.jl. It is sufficient to load your plotting package, i.e.
  Makie.jl or MeshCat.jl, **before** calling `Polyhedra.Mesh`