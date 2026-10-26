# Polyhedra.jl upcoming release notes

## Deprecations
- The dimension query `Polyhedra.FullDim(p)` is now
  `Polyhedra.typed_fulldim(p)`. Implementations of the Polyhedra interface
  should extend `typed_fulldim`; `FullDim` remains the dimension type alias.
  Methods of `Polyhedra.FullDim` defined by existing implementations are
  still used by `typed_fulldim` as a fallback.

## Breaking changes
- The undefined `start*`, `done*`, and `next*` iterator names are no longer
  exported. Use Julia's `iterate` interface instead.

## New features
- `issubset(p, q)`, i.e., `p ⊆ q`, is now implemented when `p` is a
  representation or a polyhedron and `q` is an H-representation or a
  polyhedron (see [#303]).
- `eliminate(p, delset, BlockElimination())` and
  `project(p, pset, BlockElimination())` are now implemented for any library.
  The extreme rays of the cone of multipliers are computed with the library of
  `p`. As the result may contain many redundant halfspaces, the default
  elimination algorithm only uses it if the library implements it natively.

## Bug fixes
- `issubset(h::HRepresentation, el::HRepElement, solver)` no longer throws a
  `MethodError`.
- `dualfullspace` of a `LiftedHRepresentation` now contains the origin.

[#303]: https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/303

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