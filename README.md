# Polyhedra

Polyhedra provides an unified interface for Polyhedra Manipulation Libraries such as [CDDLib.jl](https://github.com/blegat/CDDLib.jl).
These manipulation notably include the transformation from (resp. to) an inequality representation of a polyhedron to (resp. from) its generator representation (convex hull of points + conic hull of rays) and projection/elimination of a variable with e.g. [Fourier-Motzkin](https://en.wikipedia.org/wiki/Fourier%E2%80%93Motzkin_elimination).

It defines the abstract type `Polyhedron` and splits the operations on this type in two categories:

* Mandatory: Operations that needs to be implemented by the Polyhedra Manipulation Libraries: e.g. Transformation between the two representations described above and variable elimination.
* Optional: Operations that can be implemented using the other operations and hence have a default implementation: e.g. 3D visualization of the polyhedron using [GLVisualize.jl](https://github.com/JuliaGL/GLVisualize.jl), intersection of polyhedra, [Minkowski addition](https://en.wikipedia.org/wiki/Minkowski_addition) of polyhedra, ...

[![Build Status](https://travis-ci.org/blegat/Polyhedra.jl.svg?branch=master)](https://travis-ci.org/blegat/Polyhedra.jl)
