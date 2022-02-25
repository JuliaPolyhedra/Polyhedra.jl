using JuMP
import HiGHS
lp_solver = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
