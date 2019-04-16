# Optimization

A polyhedron can represents the feasible set of an optimization program.
The program is infeasible when the polyhedron is empty.

```@docs
isempty
```

If the V-representation of the polyhedron has been computed, it can be used to solve the linear program.
```@docs
VRepOptimizer
```

Otherwise, any programming solver implementing the [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl) interface can be used.
See [here](http://www.juliaopt.org/JuMP.jl/dev/installation/#Getting-Solvers-1) for a list of available solvers.
```@docs
Polyhedra.default_solver
Polyhedra.linear_objective_solver
```

## Using a polyhedron for in an optimization model

A polyhedron or representation can be used in the constraint of a JuMP model.
For instance, consider the 1-simplex:
```jldoctest jump-in-hrep
julia> using Polyhedra

julia> simplex = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)
H-representation Polyhedra.Intersection{Int64,Array{Int64,1},Int64}:
1-element iterator of HyperPlane{Int64,Array{Int64,1}}:
 HyperPlane([1, 1], 1),
2-element iterator of HalfSpace{Int64,Array{Int64,1}}:
 HalfSpace([-1, 0], 0)
 HalfSpace([0, -1], 0)
```
and the following JuMP model with two variables
```jldoctest jump-in-hrep
julia> using JuMP

julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, λ[1:2])
2-element Array{VariableRef,1}:
 λ[1]
 λ[2]
```

The variables can be constrained to belong to the simplex as follows:
```jldoctest jump-in-hrep
julia> @constraint(model, λ in simplex)
[λ[1], λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))
```
but a vector of affine or quadratic expressions can also be constrained to belong to the simplex:
```jldoctest jump-in-hrep
julia> A = [1  1
            1 -1]
2×2 Array{Int64,2}:
 1   1
 1  -1

julia> @constraint(model, A * λ in simplex)
[λ[1] + λ[2], λ[1] - λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))
```
We can verify that the model contains both constraints:
```julia
julia> model
A JuMP Model
Feasibility problem with:
Variables: 2
`Array{JuMP.VariableRef,1}`-in-`Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}`: 1 constraint
`Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},1}`-in-`Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}`: 1 constraint
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.
Names registered in the model: λ
```
When the model is solved, the constraint is automatically transformed into
appropriate constraints if the optimizer does not support consraints with the
set `Polyhedra.PolyhedraOptSet`.
```julia
julia> import GLPK

julia> optimize!(model, with_optimizer(GLPK.Optimizer))

julia> termination_status(model)
OPTIMAL::TerminationStatusCode = 1

julia> value.(λ)
2-element Array{Float64,1}:
 0.5
 0.5
```
For instance, GLPK, does not support
`Polyhedra.PolyhedraOptSet` constraints but supports `MOI.EqualTo{Float64}` and
`MOI.LessThan{Float64}`. The polyhedral constraints are therefore
bridged into several `MOI.EqualTo{Float64}` and `MOI.LessThan{Float64}`
constraints using the following
[constraint bridge](http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual/#Constraint-bridges-1):
```@docs
Polyhedra.PolyhedraToLPBridge
```

See [Polyhedral Function](https://github.com/JuliaPolyhedra/Polyhedra.jl/blob/master/examples/Polyhedral%20Function.ipynb) for an example notebook.

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
