# Optimization

A polyhedron can represents the feasible set of an optimization program.
The program is infeasible when the polyhedron is empty.

```@docs
isempty
Polyhedra.MathProgBase.linprog
```

If the V-representation of the polyhedron has been computed, it can be used to solve the linear program.
```@docs
VRepSolver
```

Otherwise, any programming solver implementing the [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface can be used. See [here](http://www.juliaopt.org/) for a list of available solvers.
```@docs
Polyhedra.default_solver
Polyhedra.solver
```

## Creating a polyhedron from the feasible set of a JuMP model

A typical application of polyhedral computation is the computation of the set of extreme points and rays of the feasible set of an optimization problem.
This comes from the fact that given a minimization of a concave function (or maximization of a convex function) on a convex feasible set (e.g. Linear Programming),
we are either in the following three situations:

- The feasible set is empty, i.e. the problem is infeasible.
- An extreme ray is optimal, i.e. the problem is unbounded (or it may also be bounded if the objective is constant along the ray).
- An extreme point is optimal.

A JuMP model is treated by `polyhedron` just like any H-representation. For example, the hypercube of dimension `n` can be created as follows
```julia
m = Model()
@variable(m, 0 ≤ x[1:n] ≤ 1)

poly = polyhedron(m, CDDLib.Library(:exact))
```

In fact, the MathProgBase representation of the feasible set of a linear program:

```math
\begin{align*}
  lb \leq Ax \leq ub\\
  l \leq x \leq u\\
\end{align*}
```

has `LPHRepresentation` as a corresponding H-representation.
A JuMP model can be converted to this representation using `LPHRepresentation(m)`.
