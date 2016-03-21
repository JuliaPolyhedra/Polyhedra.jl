.. _polyhedra-installation:

------------
Installation
------------

This section show how to install Julia, Polyhedra
and a Polyhedra Manipulation Library of your choice.

Getting Julia
^^^^^^^^^^^^^

The first step is to install Julia.
At the time of writing, the latest release of Julia is version ``0.4`` adn the version ``0.5`` is in development.
Polyhedra can be used with either Julia 0.4 and 0.5.
Download links and more detailed instructions are available on the `Julia website <http://julialang.org>`_.

Getting Polyhedra
^^^^^^^^^^^^^^^^^

Julia has a package manager that makes the installation of new packages ridiculously simple.
Open a Julia console (e.g. enter ``julia`` at the command line) and write::

    julia> Pkg.add("JuMP")

To start using Polyhedra, you can now just write::

    julia> using Polyhedra

Of course without installing a library, you won't be able to do much. See the next section on installing a library.

Getting Libraries
^^^^^^^^^^^^^^^^^

.. _polyhedra-librarytable:

+-------------------------------------------------------------+----------------------------------------------------+--------------+---------+----------------+----------------+
| Solver                                                      | Julia Package                                      | Library      | License | Exact Rational | Floating point |
+=============================================================+====================================================+==============+=========+================+================+
| `cdd <https://www.inf.ethz.ch/personal/fukudak/cdd_home/>`_ | `CDDLib.jl <https://github.com/blegat/CDDLib.jl>`_ | ``CDDLib()`` |  GPL    |  X  |  X  |
+-------------------------------------------------------------+----------------------------------------------------+--------------+---------+----------------+----------------+
| `lrs <http://cgm.cs.mcgill.ca/~avis/C/lrs.html>`_           | `LRSLib.jl <https://github.com/blegat/LRSLib.jl>`_ | ``LRSLib()`` |  GPL    |  X  |  X  |
+-------------------------------------------------------------+----------------------------------------------------+--------------+---------+----------------+----------------+
