using JuMP

import GLPK
lp_solver = with_optimizer(GLPK.Optimizer)
