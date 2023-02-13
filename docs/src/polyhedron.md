# Polyhedron

As seen in the previous section, a polyhedron can be described in 2 ways: either using the H-representation (intersection of halfspaces) or the V-representation (convex hull of points and rays).
The problem of computing the H-representation from the V-representation (or vice versa) is called the *representation conversion problem*.
It can be solved by the Double-Description method
```@docs
doubledescription
```
However, other methods exist such as the reverse search implemented by [LRS](https://github.com/JuliaPolyhedra/LRSLib.jl) and the quick hull algorithm implemented by [qhull](https://github.com/JuliaPolyhedra/QHull.jl).

This motivates the creation of a type representing polyhedra, transparently handling the conversion from H-representation to V-representation when needed for some operation.
Just like the abstract type `AbstractArray{N,T}` represents an `N`-dimensional array with elements of type `T`,
the abstract type `Polyhedron{N,T}` represents an `N`-dimensional polyhedron with elements of coefficient type `T`.

There is typically one concrete subtype of `Polyhedron` by library.
For instance, the CDD library defines `CDDLib.Polyhedron` and the LRS library defines `LRSLib.Polyhedron`.
It must be said that the type `T` is not necessarily how the elements are stored internally by the library but the polyhedron will behave just like it is stored that way.
For instance, when retreiving an H-(or V-)representation, the representation will be of type `T`.
Therefore using `Int` for `T` may result in `InexactError`.
For this reason, by default, the type `T` chosen is not a subtype of `Integer`.

A polyhedron can be created from a representation and a library using the `polyhedron` function.
```@docs
polyhedron
```

To illustrate the usage of the `polyhedron` function, consider the following representations:
```julia
hr = HalfSpace([1, 1], 1) ∩ HalfSpace([1, -1], 0) ∩ HalfSpace([-1, 0], 0)
vre = convexhull([0, 0], [0, 1], [1//2, 1//2])
vrf = convexhull([0, 0], [0, 1], [1/2, 1/2])
```

One can use the CDD library, to create an instance of a concrete subtype of `Polyhedron` as follows:
```julia
julia> using CDDLib
julia> polyf = polyhedron(hr, CDDLib.Library())
julia> typeof(polyhf)
CDDLib.CDDLib.Polyhedron{2,Float64}
```

We see that the library has choosen to deal with floating point arithmetic.
This decision does not depend on the type of `hr` but only on the instance of
`CDDLib.Library` given. `CDDLib.Library` creates `CDDLib.Polyhedron` of type
either `Float64` or `Rational{BigInt}`. One can choose the first one using
`CDDLib.Library(:float)` and the second one using `CDDLib.Library(:exact)`, by
default it is `:float`.
```julia
julia> poly = polyhedron(hr, CDDLib.Library(:exact))
julia> typeof(poly)
CDDLib.Polyhedron{2,Rational{BigInt}}
```

The first polyhedron `polyf` can also be created from its V-representation using either of the 4 following lines:
```julia
julia> polyf = polyhedron(vrf, CDDLib.Library(:float))
julia> polyf = polyhedron(vrf, CDDLib.Library())
julia> polyf = polyhedron(vre,  CDDLib.Library(:float))
julia> polyf = polyhedron(vre,  CDDLib.Library())
```

and `poly` using either of those lines:
```julia
julia> poly = polyhedron(vrf, CDDLib.Library(:exact))
julia> poly = polyhedron(vre, CDDLib.Library(:exact))
```

Of course, creating a representation in floating points with exact arithmetic works here because we have `0.5` which is `0.1` in binary but in general, is not a good idea.
```julia
julia> Rational{BigInt}(1/2)
1//2
julia> Rational{BigInt}(1/3)
6004799503160661//18014398509481984
julia> Rational{BigInt}(1/5)
3602879701896397//18014398509481984
```

## Retrieving a representation

One can retrieve an H-representation (resp. V-representation) from a polyhedron using `hrep` (resp. `vrep`).
The concrete subtype of `HRepresentation` (resp. `VRepresentation`) returned is not necessarily the same that the one used to create the polyhedron.
As a rule of thumb, it is the representation the closest to the internal representation used by the library.
```julia
julia> hr = hrep(poly)
julia> typeof(hr)
Polyhedra.LiftedHRepresentation{2,Rational{BigInt}}
julia> hr = MixedMatHRep(hr)
julia> typeof(hr)
Polyhedra.MixedMatHRep{2,Rational{BigInt}}
julia> hr.A
3x2 Array{Rational{BigInt},2}:
  1//1   1//1
  1//1  -1//1
 -1//1   0//1
julia> hr.b
3-element Array{Rational{BigInt},1}:
 1//1
 0//1
 0//1
julia> vr = vrep(poly)
julia> typeof(vr)
Polyhedra.LiftedVRepresentation{2,Rational{BigInt}}
julia> vr = MixedMatVRep(vrep)
julia> typeof(vr)
Polyhedra.MixedMatVRep{2,Rational{BigInt}}
julia> vr.V
3x2 Array{Rational{BigInt},2}:
 1//2  1//2
 0//1  1//1
 0//1  0//1

julia> vr.R
0x2 Array{Rational{BigInt},2}
```

## Checking if a representation has been computed

```@docs
hrepiscomputed
vrepiscomputed
```

## Incidence

Elements can be accessed in a representation or polyhedron using indices and `Base.get`:
```@docs
Polyhedra.Index
Polyhedra.Indices
```
The list of indices can be obtained using, e.g., `eachindex(points(rep))`.
For instance, the following prints all points using indices
```julia
for pi in eachindex(points(rep))
    @show get(rep, pi)
end
```

A point ``p`` (resp. ray ``r``) is incident to an halfspace ``\langle a, x \rangle \le \beta`` if ``\langle a, p \rangle = \beta`` (resp. ``\langle a, r \rangle = \beta``).

```@docs
incidenthalfspaces
incidenthalfspaceindices
incidentpoints
incidentpointindices
incidentrays
incidentrayindices
```

In a polyhedron, all points and rays are incident to all hyperplanes and all halfspaces are incident to all lines.
The following methods are therefore redundant, e.g. `incidenthyperplanes(p, idx)` is equivalent to `hyperplanes(p)` and `incidenthyperplaneindices(p, idx)` is equivalent to `eachindex(hyperplanes(p))`.
The methods are hence only defined for consistency.

```@docs
incidenthyperplanes
incidenthyperplaneindices
incidentlines
incidentlineindices
```

## Default libraries

The following functions allows to select a default library:
```@docs
default_library
similar_library
library
default_type
```

The following libraries serves as fallback:
```@docs
DefaultLibrary
IntervalLibrary
```

The type and library of the polyhedron obtained after applying an operation of several polyhedra (of possibly different type and/or library) is determined by the `similar` function.
```@docs
similar
```
