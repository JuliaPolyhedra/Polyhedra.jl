# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

cdd = false && try_import(:CDDLib)
lrs = false && try_import(:LRSLib)

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

# Load an available solver
lp_solver = nothing
glp = try_import(:GLPKMathProgInterface)
glp && (lp_solver = GLPKMathProgInterface.GLPKSolverLP())
if lp_solver === nothing
    cbc = try_import(:Cbc)
    if cbc; import Clp; end
    cbc && (lp_solver = Clp.ClpSolver())
end
if lp_solver === nothing
    grb = try_import(:Gurobi) # Gurobi creates BigFloat when the input is BigInt and then cannot handle it
    grb && (lp_solver = Gurobi.GurobiSolver(OutputFlag=0))
end
if lp_solver === nothing
    cpx = try_import(:CPLEX)
end
if lp_solver === nothing
    xpr = try_import(:Xpress)
    xpr && (lp_solver = Xpress.XpressSolver(OUTPUTLOG=0))
end
if lp_solver === nothing
    mos = try_import(:Mosek)
    mos && (lp_solver = Mosek.MosekSolver(LOG=0))
end
if lp_solver === nothing
    ipt = try_import(:Ipopt)
    ipt && (lp_solver = Ipopt.IpoptSolver(print_level=0))
end
if lp_solver === nothing
    eco = try_import(:ECOS)
    eco && (lp_solver = ECOS.ECOSSolver(verbose=false))
end
if lp_solver === nothing
    scs = try_import(:SCS)
    scs && (lp_solver = SCS.SCSSolver(eps=1e-6,verbose=0))
    @assert lp_solver !== nothing
end
lpsolver = tuple(lp_solver)
