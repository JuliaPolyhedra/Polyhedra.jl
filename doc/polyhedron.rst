.. _polyhedra-polyhedron:

----------
Polyhedron
----------

As seen in the previous section, a polyhedron can be described in 2 ways: either using the H-representation or the V-representation.
A typical problem is: Given the H-(or V-)representation of one or several polyhedra, what is the H-(or V-)representation of some polyhedra obtained after some operations of those initial polyhedra.
This description is similar to the description usually given to algorithms except that in that case we talk about numbers given in their binary representation and not polyhedra given in their H-(or V-)representation.
This motivates the creation of a type representing polyhedra.
Just like the abstract type ``AbstractArray{N,T}`` represents an ``N``-dimensional array of elements of type ``T``,
the abstract type ``Polyhedron{N,T}`` represents an ``N``-dimensional polyhedron of elements of type ``T``.

There is typically one concrete subtype of ``Polyhedron`` by library.
For instance, the CDD library defines ``CDDPolyhedron`` and the LRS library defines ``LRSPolyhedron``.
It must be said that the type ``T`` is not necessarily how the elements are stored internally by the library but the polyhedron will behave just like it is stored that way.
For instance, when retreiving an H-(or V-)representation, the representation will be of type ``T``.
Therefore ``Int`` for ``T`` is may result in ``InexactError``.
For this reason, by default, the type ``T`` chosen is not a subtype of ``Integer``.

Consider the representations ``hrep``, ``vrep`` and ``vrepf`` created in the preceding section.
One can use the CDD library, to create an instance of a concrete subtype of ``Polyhedron``::

    julia> polyf = polyhedron(hrep, CDDLibrary())
    julia> typeof(polyhf)
    CDDLib.CDDPolyhedron{2,Float64}

We see that the library has choosen to deal with floating point arithmetic.
This decision does not depend on the type of ``hrep`` but only on the instance of ``CDDLibrary`` given.
``CDDLibrary`` creates ``CDDPolyhedron`` of type either ``Float64`` or ``Rational{BigInt}``.
One can choose the first one using ``CDDLibrary(:float)`` and the second one using ``CDDLibrary(:exact)``, by default it is ``:float``.::

    julia> poly = polyhedron(hrep, CDDLibrary(:exact))
    julia> typeof(poly)
    CDDLib.CDDPolyhedron{2,Rational{BigInt}}


The first polyhedron ``polyf`` can also be created from its V-representation using either of the 4 following lines::

    julia> polyf = polyhedron(vrepf, CDDLibrary(:float))
    julia> polyf = polyhedron(vrepf, CDDLibrary())
    julia> polyf = polyhedron(vrep,  CDDLibrary(:float))
    julia> polyf = polyhedron(vrep,  CDDLibrary())

and ``poly`` using either of those lines::

    julia> poly = polyhedron(vrepf, CDDLibrary(:exact))
    julia> poly = polyhedron(vrep , CDDLibrary(:exact))

of course, creating a representation in floating points with exact arithmetic works here because we have ``0.5`` which is ``0.1`` in binary but in general, is not a good idea.::

    julia> Rational{BigInt}(1/2)
    1//2
    julia> Rational{BigInt}(1/3)
    6004799503160661//18014398509481984
    julia> Rational{BigInt}(1/5)
    3602879701896397//18014398509481984

Retrieving a representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can retrieve an H-representation (resp. V-representation) from a polyhedron using ``hrep`` (resp. ``vrep``).
The concrete subtype of ``HRepresentation`` (resp. ``VRepresentation``) returned is not necessarily the same that the one used to create the polyhedron.
As a rule of thumb, it is the representation the closest to the internal representation used by the library.::

    julia> hrep = hrep(poly)
    julia> typeof(hrep)
    Polyhedra.LiftedHRepresentation{2,Rational{BigInt}}
    julia> hrep = SimpleHRepresentation(hrep)
    julia> typeof(hrep)
    Polyhedra.SimpleHRepresentation{2,Rational{BigInt}}
    julia> hrep.A
    3x2 Array{Rational{BigInt},2}:
      1//1   1//1
      1//1  -1//1
     -1//1   0//1
    julia> hrep.b
    3-element Array{Rational{BigInt},1}:
     1//1
     0//1
     0//1
    julia> vrep = vrep(poly)
    julia> typeof(vrep)
    Polyhedra.LiftedVRepresentation{2,Rational{BigInt}}
    julia> vrep = SimpleVRepresentation(vrep)
    julia> typeof(vrep)
    Polyhedra.SimpleVRepresentation{2,Rational{BigInt}}
    julia> vrep.V
    3x2 Array{Rational{BigInt},2}:
     1//2  1//2
     0//1  1//1
     0//1  0//1

    julia> vrep.R
    0x2 Array{Rational{BigInt},2}
