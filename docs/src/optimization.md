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

Otherwise, any programming solver implementing the [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) interface can be used.
See [here](http://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers-1) for a list of available solvers.
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
H-representation Polyhedra.Intersection{Int64, Vector{Int64}, Int64}:
1-element iterator of HyperPlane{Int64, Vector{Int64}}:
 HyperPlane([1, 1], 1),
2-element iterator of HalfSpace{Int64, Vector{Int64}}:
 HalfSpace([-1, 0], 0)
 HalfSpace([0, -1], 0)
```
and the following JuMP model with two variables
```jldoctest jump-in-hrep
julia> using JuMP

julia> model = Model()
A JuMP Model
├ solver: none
├ objective_sense: FEASIBILITY_SENSE
├ num_variables: 0
├ num_constraints: 0
└ Names registered in the model: none

julia> @variable(model, λ[1:2])
2-element Vector{VariableRef}:
 λ[1]
 λ[2]
```

The variables can be constrained to belong to the simplex as follows:
```jldoctest jump-in-hrep
julia> @constraint(model, λ in simplex)
[λ[1], λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64, Polyhedra.Intersection{Int64, Vector{Int64}, Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))
```
but a vector of affine or quadratic expressions can also be constrained to belong to the simplex:
```jldoctest jump-in-hrep
julia> A = [1  1
            1 -1]
2×2 Matrix{Int64}:
 1   1
 1  -1

julia> @constraint(model, A * λ in simplex)
[λ[1] + λ[2], λ[1] - λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64, Polyhedra.Intersection{Int64, Vector{Int64}, Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))
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

julia> set_optimizer(model, GLPK.Optimizer)

julia> optimize!(model)

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
[constraint bridge](http://jump.dev/MathOptInterface.jl/stable/apimanual/#Constraint-bridges-1):
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
```@example lphrep
using JuMP, Polyhedra
model = Model()
@variable(model, 0 ≤ x[1:2] ≤ 1)
h = hrep(model)
```
The name of the variables for each dimension can be recovered as follows
```@example lphrep
dimension_names(h)
```
Note that the names of the JuMP variables are lost in the conversion to a
polyhedron that does not support names, e.g.,
```julia
poly = polyhedron(model, CDDLib.Library(:exact))
```
However, the ordering of the dimension of the polyhedron is guaranteed to
correspond to the order of the JuMP variables as listed by `all_variables`:
```@example lphrep
all_variables(model)
```
So `all_variables(model)[i]` provides the JuMP variable corresponding to the `i`th
dimension.
The reverse mapping can be constructed as follows:
```@example lphrep
Dict(v => i for (i, v) in enumerate(all_variables(model)))
```

## Using a Polyhedra Optimizer with MathOptInterface

Polyhedra Optimizers by dafault support only a few constraint types in MathOptInterface (MOI).
Apply `MOI.Bridges.full_bridge_optimizer` to a Polyhedra Optimizer to enable a broader set of constraint types, such as `VectorAffineFunction`:
see the [list](http://www.juliaopt.org/MathOptInterface.jl/dev/apimanual/#Constraints-by-function-set-pairs-1) from MOI.

As an example, consider the linear program:
```math
\[
\max\ c x \quad \text{s.t.}\ A x \leq b,\ x \geq 0
\]
```
where
```julia
A = [1 2; 3 1]
b = [1, 2]
c = [1, 1]
```

Let us solve this program with `CDDLib.Optimizer` in exact arithmetic.
To set up:

```julia
using CDDLib
using MathOptInterface
const MOI = MathOptInterface

m, n = size(A)
T = Rational{BigInt}

# Enable `VectorAffineTerm`, `VectorOfVariables`, `Nonnegatives`, `Nonpositives`
optimizer = MOI.Bridges.full_bridge_optimizer(CDDLib.Optimizer{T}(), T)

# max cx
x = MOI.add_variables(optimizer, n)
MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
        MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}.(c, x), 0))
MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)

# Ax - b <= 0
terms = MOI.VectorAffineTerm{T}.(
    1:m, MOI.ScalarAffineTerm{T}.(A, reshape(x, 1, n))
)
f = MOI.VectorAffineFunction{T}(vec(terms), -b)
MOI.add_constraint(optimizer, f, MOI.Nonpositives(m))

# x >= 0
MOI.add_constraint(optimizer, MOI.VectorOfVariables(x), MOI.Nonnegatives(n))
```

To solve:

```julia
MOI.optimize!(optimizer)
MOI.get(optimizer, MOI.VariablePrimal(), x)
```

```
2-element Array{Rational{BigInt},1}:
 3//5
 1//5
```
