## Projection/Elimination

Consider the polyhedron created in the beginning of this section. As a reminder, it represents the following H-representation:
```math
\begin{align*}
  x_1 + x_2 &\leq 1 \\
  x_1 - x_2 &\leq 0 \\
  x_1 & \geq 0.
\end{align*}
```

One can verify that for any ``0 \leq x_2 \leq 1``, there exists a value ``x_1`` such that ``(x_1, x_2)`` is in this polyhedron.
This means that the H-representation obtained by eliminating ``x_1`` is:

```math
\begin{align*}
  x_1 & \leq 1 \\
  x_1 & \geq 0.
\end{align*}
```

where ``x_1`` in the H-representation above represents ``x_2`` in the previous one.
This can be obtained as follows
```julia
julia> poly_x2 = eliminate(poly, [1])
julia> hrep(poly_x2)
H-representation
begin
 2 2 rational
 1//1 -1//1
 0//1 1//1
end
```

There is two methods of computing the elimination implemented in CDDLib: Fourier-Motzkin elimination and block elimination.
As written by K. Fukuda in CDD's documentation, "[Block elimination] might be a faster way to eliminate variables than the repeated [Fourier-Motzkin elimination] when the number of variables to eliminate is large".
You can specify the method to use as a third argument, e.g. `eliminate(poly, [1], FourierMotzkin())`, `eliminate(poly, [1], BlockElimination())`.
A third method can be chosen: `ProjectGenerators`.
It computes the V-representation and then project each of its elements.
This is the method of choice when the V-representation is already computed.
```@docs
FourierMotzkin
BlockElimination
ProjectGenerators
```

If nothing is specified as in the block of code above, the behavior depends on the polyhedral library.
If neither Fourier-Motzkin nor block elimination is implemented or if the V-representation is already computed then `:ProjectGenerators` is chosen.
Otherwise, Polyhedra lets the library decide. In CDDLib, `:FourierMotzkin` is chosen when only the last dimension needs to be eliminated and `:BlockElimination` is chosen otherwise.
Note that CDDLib only supports projecting the last trailing dimensions.

```@docs
eliminate
project
fixandeliminate
```
