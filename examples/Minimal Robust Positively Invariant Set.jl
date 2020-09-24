# # Minimal Robust Positively Invariant Set

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Minimal Robust Positively Invariant Set.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Minimal Robust Positively Invariant Set.ipynb)

# ### Introduction

# In this example, we implement the method presented in [RKKM05] to compute a robust positively invariant polytope of a linear system under a disturbance bounded by a polytopic set.
#
# We consider the example given in equation (15) of [RKKM05]:
# ```math
# x^+ =
# \begin{bmatrix}
#   1 & 1\\
#   0 & 1
# \end{bmatrix}x +
# \begin{bmatrix}
#   1\\
#   1
# \end{bmatrix} u
# + w
# ```
# with the state feedback control $u(x) = -\begin{bmatrix}1.17 & 1.03\end{bmatrix} x$.
# The controlled system is therefore
# ```math
# x^+ =
# \left(\begin{bmatrix}
#   1 & 1\\
#   0 & 1
# \end{bmatrix} -
# \begin{bmatrix}
#   1\\
#   1
# \end{bmatrix}
# \begin{bmatrix}1.17 & 1.03\end{bmatrix}\right)x
# + w =
# \begin{bmatrix}
#   -0.17 & -0.03\\
#   -1.17 & -0.03
# \end{bmatrix}x
# + w.
# ```
#
# [RKKM05] Sasa V. Rakovic, Eric C. Kerrigan, Konstantinos I. Kouramas, David Q. Mayne *Invariant approximations of the minimal robust positively Invariant set*. IEEE Transactions on Automatic Control 50 (**2005**): 406-410.

A = [1 1; 0 1] - [1; 1] * [1.17 1.03]

# The set of disturbance is the unit ball of the infinity norm.

using Test #jl
using Polyhedra
Wv = vrep([[x, y] for x in [-1.0, 1.0] for y in [-1.0, 1.0]])

# We will use the default library for this example but feel free to pick any other library from [this list of available libraries](https://juliapolyhedra.github.io/) such as [CDDLib](https://github.com/JuliaPolyhedra/CDDLib.jl).
# The LP solver used to detect redundant points in the V-representation is [GLPK](https://github.com/JuliaOpt/GLPK.jl). Again, you can replace it with any other solver listed [here](http://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers-1) that supports LP.

using GLPK
using JuMP
lib = DefaultLibrary{Float64}(GLPK.Optimizer)
W = polyhedron(Wv, lib)

# The $F_s$ function of equation (2) of [RKKM05] can be implemented as follows.

function Fs(s::Integer, verbose=1)
    @assert s ≥ 1
    F = W
    A_W = W
    for i in 1:(s-1)
        A_W = A * A_W
        F += A_W
        if verbose ≥ 1
            println("Number of points after adding A^$i * W: ", npoints(F))
        end
        removevredundancy!(F)
        if verbose ≥ 1
            println("Number of points after removing redundant ones: ", npoints(F))
        end
    end
    return F
end

# We can see below that only the V-representation is computed. In fact, no H-representation was ever computed during `Fs`. Computing $AW$ is done by multiplying all the points by $A$ and doing the Minkowski sum is done by summing each pair of points. The redundancy removal is carried out by CDD's internal LP solver.

@time Fs(4) #!jl
@test npoints(Fs(4)) == 16 #jl

# The Figure 1 of [RKKM05] can be reproduced as follows:

using Plots          #!jl
plot()               #!jl
for i in 10:-1:1     #!jl
    plot!(Fs(i, 0))  #!jl
end                  #!jl
# The cell needs to return the plot for it to be displayed
# but the `for` loop returns `nothing` so we add this dummy `plot!` that returns the plot
plot!()              #!jl

# Now, suppose we want to compute an invariant set by scaling $F_s$ by the appropriate $\alpha$.
# In equation (11) of [RKKM05], we want to check whether $A^s W \subseteq \alpha W$ which is equivalent to $W \subseteq \alpha A^{-s} W$.
# Note that `A^s \ W` triggers the computation of the H-representation of `W` and `A_W` is H-represented.

function αo(s)
    A_W = A^s \ W
    hashyperplanes(A_W) && error("HyperPlanes not supported")
    return maximum([Polyhedra.support_function(h.a, W) / h.β for h in halfspaces(A_W)])
end

# We obtain $\alpha \approx 1.9 \cdot 10^{-5}$ like in [RKKM05].

α = αo(10)
@test α ≈ 1.91907e-5 #jl

# The scaled set is is the following:

using Plots                    #!jl
plot((1 - α)^(-1) * Fs(10, 0)) #!jl
