# Inspired from MathOptInterface/src/Test/config.jl
"""
    @polytestset setname subsets

Defines a function `setnametest(lib, exclude)` that runs the tests defined in the dictionary `setnametests`
with the library `lib` except the tests whose dictionary key is in `exclude`.
If `subsets` is `true` then each test runs in fact multiple tests hence the `exclude` argument is passed
as it can also contains test to be excluded from these subsets of tests.
"""
macro polytestset(setname, subsets=false)
    testname = Symbol(string(setname) * "test")
    testdict = Symbol(string(testname) * "s")
    if subsets
        runtest = :( f(lib, exclude) )
    else
        runtest = :( f(lib) )
    end
    esc(:(
        function $testname(lib::Polyhedra.Library, exclude::Vector{String} = String[])
            for (name,f) in $testdict
                if name in exclude
                    continue
                end
                @testset "$name" begin
                    $runtest
                end
            end
        end
    ))
end
