import GLPK
# Need `"presolve" => GLPK.ON` for `detect_new_linearities`, see
# https://travis-ci.org/github/JuliaPolyhedra/Polyhedra.jl/jobs/691916637#L486
lp_solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)
