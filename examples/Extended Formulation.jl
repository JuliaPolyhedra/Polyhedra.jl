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

# In the following, `diagonal` and `antidiag` are projections of extended Hrepresentations of dimension 7
# hence `diamond` is a projection of an extended H-representation of dimensions 12.

diagonal = convexhull(translate(square, [-2, -2]), translate(square, [2, 2]))
@test fulldim(diagonal.set) == 7 #jl
antidiag = convexhull(translate(square, [-2,  2]), translate(square, [2, -2]))
@test fulldim(antidiag.set) == 7 #jl
diamond = diagonal ∩ antidiag
@test fulldim(diamond.set) == 12 #jl

# We don't need to compute the result of the projection to solve an optimization problem
# over `diamond`.

import GLPK
model = Model(GLPK.Optimizer)
@variable(model, x[1:2])
@constraint(model, x in diamond)
@objective(model, Max, x[2])
optimize!(model)
value.(x) #!jl
@test value.(x) ≈ [0.0, 2.0] #jl

# In the optimization problem, above, the auxiliary variables of the extended formulation
# are transparently added inside a bridge.
# To manipulate the auxiliary variables, one can use the extended H-representation directly instead of its projection.

import GLPK
model = Model(GLPK.Optimizer)
@variable(model, x[1:12])
@constraint(model, x in diamond.set)
@objective(model, Max, x[2])
optimize!(model)
value.(x) #!jl
@test value.(x) ≈ [0.0, 2.0, -0.75, -0.25, 0.75, 2.25, 0.25, -0.75, 2.25, 0.75, -0.25, 0.75] #jl
