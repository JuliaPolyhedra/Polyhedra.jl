# # Extended formulation

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Extended Formulation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Extended Formulation.ipynb)

# In this notebook, we show how to work with the extended formulation of a polyhedron.
# The convex hull of the union of polyhedra that are H-represented can be obtained
# as the projection of a H-representation [Theorem 3.3, B85].
# In order to use the resulting polyhedron as a constraint set in an optimization problem,
# there is no need to compute the resulting H-representation of this projection.
# Moreover, other operations such as intersection are also implemented between extended H-representations.
# We illustrate this with a simple example. We start by defining the H-representation of a square with JuMP.

# [B85] Balas, E., 1985.
# *Disjunctive programming and a hierarchy of relaxations for discrete optimization problems*.
# SIAM Journal on Algebraic Discrete Methods, 6(3), pp.466-486.

using JuMP
model = Model()
@variable(model, -1 <= x <= 1)
@variable(model, -1 <= y <= 1)
using Polyhedra
square = hrep(model)

# Note that the names of the JuMP variables are used as the names of the corresponding dimensions for
# the polyhedron `square`.

dimension_names(square) #!jl
@test dimension_names(square) == ["x", "y"] #jl

# In the following, `diagonal` and `antidiag` are projections of extended H-representations of dimension 7
# hence `diamond` is a projection of an extended H-representation of dimensions 12.

diagonal = convexhull(translate(square, [-2, -2]), translate(square, [2, 2]))
@test fulldim(diagonal.set) == 7 #jl
@test dimension_names(diagonal.set) == ["x"; "y"; fill("", 5)] #jl
antidiag = convexhull(translate(square, [-2,  2]), translate(square, [2, -2]))
@test fulldim(antidiag.set) == 7 #jl
@test dimension_names(antidiag.set) == ["x"; "y"; fill("", 5)] #jl
diamond = diagonal ∩ antidiag
@test fulldim(diamond.set) == 12 #jl
@test dimension_names(diamond.set) == ["x"; "y"; fill("", 10)] #jl

# Note that the names the first two dimensions are still identical to the names of the JuMP variables
# and the auxiliary variables have no name.

dimension_names(diamond.set) #!jl
@test dimension_names(diamond.set) == ["x"; "y"; fill("", 10)] #jl

# We don't need to compute the result of the projection to solve an optimization problem
# over `diamond`. For instance, to compute the maximal value that `y` can take over this
# polytope with the GLPK solver, we can do as follows.
# Note that if we use anonymous JuMP variables, the name of the JuMP variables will be
# the names of the corresponding dimensions of the polyhedron.
# Therefore, we can retrieve the JuMP variable according to the corresponding dimension name
# with `variable_by_name`.

import GLPK
model = Model(GLPK.Optimizer)
@variable(model, [1:2] in diamond)
x = variable_by_name(model, "x")
y = variable_by_name(model, "y")
@objective(model, Max, y)
optimize!(model)
value(x), value(y) #!jl
@test value(x) ≈ 0.0 atol=1e-6 #jl
@test value(y) ≈ 2.0           #jl

# In the optimization problem, above, the auxiliary variables of the extended formulation
# are transparently added inside a bridge.
# To manipulate the auxiliary variables, one can use the extended H-representation directly instead of its projection.
# Note that as the auxiliary dimensions have no name, we cannot use `variable_by_name` to retrieve the corresponding
# JuMP variables. We can instead catch the returned value of `@variable` in some variable `v` in order to use anonymous JuMP
# variables while still assigning the created JuMP variables to `v`.

import GLPK
model = Model(GLPK.Optimizer)
v = @variable(model, [1:12] in diamond.set)
y = variable_by_name(model, "y")
@objective(model, Max, y)
optimize!(model)
value.(v) #!jl
@test value.(v) ≈ [0.0, 2.0, -0.75, -0.25, 0.75, 2.25, 0.25, -0.75, 2.25, 0.75, -0.25, 0.75] #jl
