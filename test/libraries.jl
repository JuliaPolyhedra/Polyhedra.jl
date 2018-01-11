# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

cdd = try_import(:CDDLib)
lrs = try_import(:LRSLib)

# Create library lists
libraries = PolyhedraLibrary[]
push!(libraries, SimplePolyhedraLibrary{Rational{BigInt}}())
push!(libraries, SimplePolyhedraLibrary{Float64}())
cdd && push!(libraries, CDDLib.CDDLibrary(:float))
cdd && push!(libraries, CDDLib.CDDLibrary(:exact))
lrs && push!(libraries, LRSLib.LRSLibrary())
# Floating point arithmetic libraries
float_libraries = PolyhedraLibrary[]
push!(float_libraries, SimplePolyhedraLibrary{Float64}())
cdd && push!(float_libraries, CDDLib.CDDLibrary(:float))
# Exact arithmetic libraries
exact_libraries = PolyhedraLibrary[]
push!(exact_libraries, SimplePolyhedraLibrary{Rational{BigInt}}())
cdd && push!(exact_libraries, CDDLib.CDDLibrary(:exact))
lrs && push!(exact_libraries, LRSLib.LRSLibrary())

# Load available solvers
glp = try_import(:GLPKMathProgInterface)
cbc = try_import(:Cbc)
if cbc; import Clp; end
grb = false && try_import(:Gurobi) # Gurobi creates BigFloat when the input is BigInt and then cannot handle it
cpx = try_import(:CPLEX)
xpr = try_import(:Xpress)
mos = try_import(:Mosek)
ipt = try_import(:Ipopt)
eco = try_import(:ECOS)
scs = try_import(:SCS)

# Create solver lists
lp_solvers = Any[]
grb && push!(lp_solvers, Gurobi.GurobiSolver(OutputFlag=0))
xpr && push!(lp_solvers, Xpress.XpressSolver(OUTPUTLOG=0))
mos && push!(lp_solvers, Mosek.MosekSolver(LOG=0))
cbc && push!(lp_solvers, Clp.ClpSolver())
glp && push!(lp_solvers, GLPKMathProgInterface.GLPKSolverLP())
ipt && push!(lp_solvers, Ipopt.IpoptSolver(print_level=0))
eco && push!(lp_solvers, ECOS.ECOSSolver(verbose=false))
scs && push!(lp_solvers, SCS.SCSSolver(eps=1e-6,verbose=0))
