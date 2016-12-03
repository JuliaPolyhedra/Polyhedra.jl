.. _polyhedra-installation:

------------
Installation
------------

This section shows how to install Julia, Polyhedra
and a Polyhedra Manipulation Library of your choice.

Getting Julia
^^^^^^^^^^^^^

The first step is to install Julia.
At the time of writing, the latest release of Julia is version ``0.4`` and the version ``0.5`` is in development.
Polyhedra can be used with either Julia ``0.4`` and ``0.5``.
Download links and more detailed instructions are available on the `Julia website <http://julialang.org>`_.

Getting Polyhedra
^^^^^^^^^^^^^^^^^

Julia has a package manager that makes the installation of new packages ridiculously simple.
Open a Julia console (e.g. enter ``julia`` at the command line) and write::

    julia> Pkg.add("Polyhedra")

To start using Polyhedra, you can now just write::

    julia> using Polyhedra

Of course without installing a library, you won't be able to do much. See the next section on installing a library.

Getting Libraries
^^^^^^^^^^^^^^^^^

.. _polyhedra-librarytable:

Many C libraries are are available for manipulating Polyhedra.
Some of them works with floating point arithmetic and some of them can do the computation exactly using rational arithmetic and multiple precision libraries such as `GMP <https://gmplib.org/>`_.
Julia also natively support Rational arithmetic using multiple precision libraries and of course floating point arithmetic.
That makes the use of both arithmetic very easy and transparent.

The following table provides a list of Polyhedra Manipulation Libraries.
When they have a Julia library implementing the interface of ``Polyhedra.jl`` then the "Library" column shows the name of the library.

+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| Solver                                                               | Julia Package                                                   | Library             | License | Exact Rational | Floating point |
+======================================================================+=================================================================+=====================+=========+================+================+
| `cdd <https://www.inf.ethz.ch/personal/fukudak/cdd_home/>`_          | `CDDLib.jl <https://github.com/blegat/CDDLib.jl>`_              | ``CDDLib()``        |  GPL    |        X       |        X       |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `lrs <http://cgm.cs.mcgill.ca/~avis/C/lrs.html>`_                    | `LRSLib.jl <https://github.com/blegat/LRSLib.jl>`_              | ``LRSLib()``        |  GPL    |        X       |        X       |
+----------------------------------------------------------------------+-----------------------------------------------------------------+--------------+-------------------+-------------+----------------+
| `ConvexHull <https://github.com/joehuchette/ConvexHull.jl>`_         | `ConvexHull.jl <https://github.com/joehuchette/ConvexHull.jl>`_ | ``ConvexHullLib()`` |  MIT    |        X       |                |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `qhull <http://www.qhull.org/>`_                                     | `QHull.jl <https://github.com/davidavdav/QHull.jl>`_            |                     |         |                |        X       |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `CHull2d <https://github.com/cc7768/CHull2d.jl>`_                    | `CHull2d.jl <https://github.com/cc7768/CHull2d.jl>`_            |                     |  MIT    |        X       |        X       |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `NewPolka <http://pop-art.inrialpes.fr/people/bjeannet/newpolka/>`_  | None                                                            |                     |  GPL    | X              |                |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `pd <http://www.cs.unb.ca/~bremner/pd/>`_                            | None                                                            |                     |  GPL    |        X       |                |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `Parma Polyhedra Library <http://bugseng.com/products/ppl/>`_        | None                                                            |                     |  GPL    | X              |                |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+
| `porta <http://comopt.ifi.uni-heidelberg.de/software/PORTA/>`_       | None                                                            |                     |  GPL    | X (overflow !) |                |
+----------------------------------------------------------------------+-----------------------------------------------------------------+---------------------+---------+----------------+----------------+

Please let me know if you plan to write a new wrapper (or an implementation in pure Julia).
Since libraries use different algorithms, no library is better for every problem; `here <http://cgm.cs.mcgill.ca/~avis/doc/avis/ABS96a.ps>`_ and `here <http://bugseng.com/products/ppl/performance>`_ are comparisons.
