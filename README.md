# Polyhedra

| **Documentation** | **PackageEvaluator** | **Build Status** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

[<img src="examples/drakeperm.png" height="240">](https://github.com/JuliaPolyhedra/Polyhedra.jl/tree/master/examples/drakeperm.jl)
[<img src="examples/glvizperm.png" height="240">](https://github.com/JuliaPolyhedra/Polyhedra.jl/tree/master/examples/glvizperm.jl)

Polyhedra provides an unified interface for Polyhedral Computation Libraries such as [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl).
These manipulation notably include the transformation from (resp. to) an inequality representation of a polyhedron to (resp. from) its generator representation (convex hull of points + conic hull of rays) and projection/elimination of a variable with e.g. [Fourier-Motzkin](https://en.wikipedia.org/wiki/Fourier%E2%80%93Motzkin_elimination).

Documentation can be found [here](http://polyhedra.readthedocs.org/).

It defines the abstract type `Polyhedron` and splits the operations on this type in two categories:

* Mandatory: Operations that needs to be implemented by the Polyhedral Computation Libraries: e.g. Transformation between the two representations described above and variable elimination.
* Optional: Operations that can be implemented using the other operations and hence have a default implementation: e.g. creation of the polyhedron from the feasible set of a [JuMP](https://github.com/JuliaOpt/JuMP.jl) model, linear transformation, intersection, [Minkowski addition](https://en.wikipedia.org/wiki/Minkowski_addition), decomposition into points and faces for e.g. 3D visualization using [DrakeVisualizer](https://github.com/rdeits/DrakeVisualizer.jl) or [GLVisualize.jl](https://github.com/JuliaGL/GLVisualize.jl)...

Here is a list of the Polyhedral Computation Libraries implementing the interface of `Polyhedra.jl`

* [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl): [CDD](https://www.inf.ethz.ch/personal/fukudak/cdd_home/) wrapper.
* [ConvexHull.jl](https://github.com/JuliaPolyhedra/ConvexHull.jl): Pure Julia implementation of the double description algorithm.
* [LRSLib.jl](https://github.com/JuliaPolyhedra/LRSLib.jl): [LRS](http://cgm.cs.mcgill.ca/~avis/C/lrs.html) wrapper.
* [QHull.jl](https://github.com/davidavdav/QHull.jl): [qhull](http://www.qhull.org/) wrapper.

Please let me know if you plan to write a new wrapper (or an implementation in pure Julia).
Since libraries (`cdd`, `lrs`, `pd`, `porta`, `qhull`, ...) use different algorithms, there is no better library for every problem; [here](http://cgm.cs.mcgill.ca/~avis/doc/avis/ABS96a.ps) is a comparison.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: http://polyhedra.readthedocs.io/en/stable
[docs-latest-url]: http://polyhedra.readthedocs.io/en/latest

[pkg-0.5-img]: http://pkg.julialang.org/badges/Polyhedra_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=Polyhedra
[pkg-0.6-img]: http://pkg.julialang.org/badges/Polyhedra_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=Polyhedra

[build-img]: https://travis-ci.org/JuliaPolyhedra/Polyhedra.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaPolyhedra/Polyhedra.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/q8shfbgt5jcv2g20/branch/master?svg=true
[winbuild-url]: https://ci.appveyor.com/project/blegat/polyhedra-jl/branch/master
[coveralls-img]: https://coveralls.io/repos/github/JuliaPolyhedra/Polyhedra.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaPolyhedra/Polyhedra.jl?branch=master
[codecov-img]: http://codecov.io/github/JuliaPolyhedra/Polyhedra.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaPolyhedra/Polyhedra.jl?branch=master
