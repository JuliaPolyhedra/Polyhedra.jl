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

    julia> using CDDLib
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

Creating a polyhedron from the feasible set of a JuMP model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A typical application of polyhedral computation is the computation of the set of extreme points and rays of the feasible set of an optimization problem.
This comes from the fact that given a minimization of a concave function (or maximization of a convex function) on a convex feasible set (e.g. Linear Programming),
we are either in the following three situations:

- The feasible set is empty, i.e. the problem is infeasible.
- An extreme ray is optimal, i.e. the problem is unbounded (or it may also be bounded if the objective is constant along the ray).
- An extreme point is optimal.

A JuMP model is treated by ``polyhedron`` just like any H-representation. For example, the hypercube of dimension ``n`` can be created as follows::

    m = Model()
    @variable(m, 0 ≤ x[1:n] ≤ 1)

    poly = polyhedron(m, CDDLibrary(:exact))

In fact, the MathProgBase representation of the feasible set of a linear program:

.. math::
    lb \leq Ax \leq ub\\
    l \leq x \leq u\\

has ``LPHRepresentation`` as a corresponding H-representation.
A JuMP Model can be converted to this representation using ``LPHRepresentation(m)``.

Projecting a polyhedron
^^^^^^^^^^^^^^^^^^^^^^^

Consider the polyhedron created in the beginning of this section. As a reminder, it represents the following H-representation:

.. math::

   x_1 + x_2 &\leq 1 \\
   x_1 - x_2 &\leq 0 \\
   x_1 & \geq 0.

One can verify that for any :math:`0 \leq x_2 \leq 1`, there exists a value :math:`x_1` such that :math:`(x_1, x_2)` is in this polyhedron.
This means that the H-representation obtained by eliminating :math:`x_1` is:

.. math::

   x_1 & \leq 1 \\
   x_1 & \geq 0.

where :math:`x_1` in the H-representation above represents :math:`x_2` in the previous one.
This can be obtained as follows:

    julia> poly_x2 = eliminate(poly, IntSet([1]))
    julia> hrep(poly_x2)
    H-representation
    begin
     2 2 rational
     1//1 -1//1
     0//1 1//1
    end

There is two methods of computing the elimination implemented in `CDDLib`: Fourier-Motzkin elimination and block elimination.
As written by K. Fukuda in CDD's documentation, "[Block elimination] might be a faster way to eliminate variables than the repeated [Fourier-Motzkin elimination] when the number of variables to eliminate is large".
You can specify the method to use as a third argument, e.g. `eliminate(poly, IntSet([1]), :FourierMotzkin)`, `eliminate(poly, IntSet([1]), :BlockElimination)`.
A third method can be chosen: `:ProjectGenerators`.
It computes the V-representation and then project each of its elements.
This is the method of choice when the V-representation is already computed.

If nothing is specified as in the block of code above, the behavior depends on the polyhedral library.
If neither Fourier-Motzkin nor block elimination is implemented or if the V-representation is already computed then `:ProjectGenerators` is chosen.
Otherwise, Polyhedra lets the library decide. In CDDLib, `:FourierMotzkin` is chosen when only the last dimension needs to be eliminated and `BlockElimination` is chosen otherwise.
Note that CDDLib only supports projecting the last trailing dimensions.
