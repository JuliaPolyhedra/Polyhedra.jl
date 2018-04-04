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
lpsolver = nothing
glp = try_import(:GLPKMathProgInterface)
glp && (lpsolver = GLPKMathProgInterface.GLPKSolverLP())
if lpsolver === nothing
    cbc = try_import(:Cbc)
    if cbc; import Clp; end
    cbc && (lpsolver = Clp.ClpSolver())
end
if lpsolver === nothing
    grb = try_import(:Gurobi) # Gurobi creates BigFloat when the input is BigInt and then cannot handle it
    grb && (lpsolver = Gurobi.GurobiSolver(OutputFlag=0))
end
if lpsolver === nothing
    cpx = try_import(:CPLEX)
end
if lpsolver === nothing
    xpr = try_import(:Xpress)
    xpr && (lpsolver = Xpress.XpressSolver(OUTPUTLOG=0))
end
if lpsolver === nothing
    mos = try_import(:Mosek)
    mos && (lpsolver = Mosek.MosekSolver(LOG=0))
end
if lpsolver === nothing
    ipt = try_import(:Ipopt)
    ipt && (lpsolver = Ipopt.IpoptSolver(print_level=0))
end
if lpsolver === nothing
    eco = try_import(:ECOS)
    eco && (lpsolver = ECOS.ECOSSolver(verbose=false))
end
if lpsolver === nothing
    scs = try_import(:SCS)
    scs && (lpsolver = SCS.SCSSolver(eps=1e-6,verbose=0))
    @assert lpsolver !== nothing
end
