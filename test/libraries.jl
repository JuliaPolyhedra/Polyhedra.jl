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
# Floating point arithmetic libraries
float_libraries = Any[]
cdd && push!(float_libraries, CDDLib.CDDLibrary(:float))
# Exact arithmetic libraries
exact_libraries = Any[]
cdd && push!(exact_libraries, CDDLib.CDDLibrary(:exact))
lrs && push!(exact_libraries, LRSLib.LRSLibrary())
