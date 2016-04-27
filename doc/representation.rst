.. _polyhedra-representation:

--------------
Representation
--------------

Polyhedra can be described in 2 different ways.

1. H-representation: As the intersection of finitely many halfspaces given by its facets.
2. V-representation: As the convex hull of its vertices + the conic hull of its rays where '+' is the Minkowski sum.

In `Polyhedra.jl`, those representations are given the respective abstract types `HRepresentation` and `VRepresentation`.

For instance, consider the 2-dimensional polyhedron described by the following H-representation:

.. math::

   x_1 + x_2 &\geq 1 \\
   x_1 - x_2 &\leq 0 \\
   x_1 & \geq 0.

This set of inequalities can be written in the matrix form :math:`Ax \leq b` where

.. math::

   A = \begin{pmatrix}1 & 1\\1 & -1\\1 & 0\end{pmatrix}, b = \begin{pmatrix}1\\0\\0\end{pmatrix}.

Let's create this H-representation using the concrete subtype `SimpleHRepresentation` of the abstract type `HRepresentation`.

    julia> A = [1 1;1 -1;1 0]
    julia> b = [1;0;0
    julia> hrep = SimpleHRepresentation(A, b)
    julia> typeof(hrep)
    Polyhedra.SimpleHRepresentation{2,Int64}
