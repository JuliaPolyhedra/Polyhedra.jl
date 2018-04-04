# Optimization

A polyhedron can represents the feasible set of an optimization program.
The program is infeasible when the polyhedron is empty.

```@docs
isempty
MathProgBase.linprog
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

The feasible set of a JuMP model can be obtained using `vrep`.
```@docs
vrep(::JuMP.Model)
```
