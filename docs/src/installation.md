# Installation

This section shows how to install Julia, Polyhedra
and a Polyhedra Manipulation Library of your choice.

## Getting Julia

The first step is to install Julia.
Polyhedra supports Julia v0.4 to Julia v0.6 but the latest version only supports Julia v0.6 (the command `Pkg.add` will automatically install the latest compatible version of Polyhedra).
Download links and more detailed instructions are available on the [Julia website](http://julialang.org).

## Getting Polyhedra

Open a Julia console (e.g. enter `julia` at the command line) and write
```julia
julia> Pkg.add("Polyhedra")
```

To start using Polyhedra, you can now just write
```julia
julia> using Polyhedra
```

Polyhedra includes a default library supporting every operation but external library can also be used.
See the next section on installing a library.

## Getting Libraries

Many C libraries are are available for manipulating Polyhedra.
Some of them works with floating point arithmetic and some of them can do the computation exactly using rational arithmetic and multiple precision libraries such as [GMP](https://gmplib.org/).
Julia also natively support Rational arithmetic using multiple precision libraries and of course floating point arithmetic.
That makes the use of both arithmetic very easy and transparent.

The following table provides a list of Polyhedra Manipulation Libraries.
When they have a Julia library implementing the interface of `Polyhedra.jl` then the "Library" column shows the name of the library.

| Solver                                                                | Julia Package                                                    | Library           | License | Exact Rational | Floating point |
|-----------------------------------------------------------------------|------------------------------------------------------------------|-------------------|---------|----------------|----------------|
| [cdd](https://www.inf.ethz.ch/personal/fukudak/cdd_home/)             | [CDDLib.jl](https://github.com/JuliaPolyhedra/CDDLib.jl)         | `CDDLibrary()`    |  GPL    |        X       |        X       |
| [ConvexHull](https://github.com/JuliaPolyhedra/ConvexHull.jl)         | [ConvexHull.jl](https://github.com/JuliaPolyhedra/ConvexHull.jl) | `ConvexHullLib()` |  MIT    |        X       |                |
| [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html)                       | [LRSLib.jl](https://github.com/JuliaPolyhedra/LRSLib.jl)         | `LRSLibrary()`    |  GPL    |        X       |        X       |
| [qhull](http://www.qhull.org/)                                        | [QHull.jl](https://github.com/davidavdav/QHull.jl)               | `QHullLib()`      |         |                |        X       |
| [CHull2d](https://github.com/cc7768/CHull2d.jl)                       | [CHull2d.jl](https://github.com/cc7768/CHull2d.jl)               |                   |  MIT    |        X       |        X       |
| [NewPolka](http://pop-art.inrialpes.fr/people/bjeannet/newpolka/)     | None                                                             |                   |  GPL    |        X       |                |
| [Parma Polyhedra Library](http://bugseng.com/products/ppl/)           | None                                                             |                   |  GPL    |        X       |                |
| [pd](http://www.cs.unb.ca/~bremner/pd/)                               | None                                                             |                   |  GPL    |        X       |                |
| [porta](http://comopt.ifi.uni-heidelberg.de/software/PORTA/)          | None                                                             |                   |  GPL    | X (overflow !) |                |

Please let me know if you plan to write a new wrapper (or an implementation in pure Julia).
Since libraries use different algorithms, no library is better for every problem; [here](http://cgm.cs.mcgill.ca/~avis/doc/avis/ABS96a.ps) and [here](http://bugseng.com/products/ppl/performance) are comparisons.
