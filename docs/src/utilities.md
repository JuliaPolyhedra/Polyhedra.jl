# Utilities

## Operations

```@docs
+
*
\
/
intersect
intersect!
convexhull
convexhull!
translate
polar
```

## Volume

```@docs
volume
surface
center_of_mass
```

## Largest inscribed ball with center

```@docs
maximum_radius_with_center
```

## Chebyshev center

```@docs
chebyshevcenter
hchebyshevcenter
vchebyshevcenter
```

## Defining new representation

The following macros make it easy to define new representations:
```@docs
Polyhedra.@subrepelem
Polyhedra.@norepelem
Polyhedra.@vecrepelem
```
