# Installation

This section shows how to install Julia, Polyhedra
and a Polyhedra Manipulation Library of your choice.

## Getting Julia

The first step is to install Julia.
Polyhedra supports Julia v0.4 to Julia v1.0 but the latest version only supports Julia v0.7 and v1.0.
Download links and more detailed instructions are available on the [Julia website](http://julialang.org).

## Getting Polyhedra

Open a Julia console (e.g. enter `julia` at the command line) and write
```julia
] add Polyhedra
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
A list of Polyhedra Manipulation Libraries is available in [the JuliaPolyhera website](https://juliapolyhedra.github.io/).
