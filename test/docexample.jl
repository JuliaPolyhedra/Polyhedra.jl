import RecipesBase
function doctest(lib::Polyhedra.Library)
    A = [1 1;1 -1;-1 0]
    b = [1,0,0]
    V = [1//2 1//2; 0 1; 0 0]
    hr = hrep(A, b)
    p = polyhedron(hr, lib)
    @test p ⊆ HalfSpace([1, -1], 0)
    @test !(p ⊆ HalfSpace([1, 1], 0))
    @test !(p ⊆ HyperPlane([-1, 0], 0))
    inequality_fulltest(p, A, b, BitSet())
    generator_fulltest(p, V)
    @test hrepiscomputed(p)
    @test vrepiscomputed(p) # computation triggered in generator_fulltest
    @test [1//4, 1//2] in p
    @test ininterior([1//4, 1//2], p)
    @test inrelativeinterior([1//4, 1//2], p)
    @test !([1, 1//2] in p)
    @test !ininterior([1, 1//2], p)
    @test !inrelativeinterior([1, 1//2], p)
    @test [0, 1] in p
    @test !ininterior([0, 1], p)
    @test !inrelativeinterior([0, 1], p)
    @test !(Ray([0, 1]) in p)
    @test !ininterior(Ray([0, 1]), p)
    @test !inrelativeinterior(Ray([0, 1]), p)
    Ap = reshape([1, -1], 2, 1)
    bp = [1, 0]
    inequality_fulltest(eliminate(p, [1]), Ap, bp, BitSet())
    inequality_fulltest(project(p, [2]), Ap, bp, BitSet())
    inequality_fulltest(eliminate(p, [1], ProjectGenerators()), Ap, bp, BitSet())
    inequality_fulltest(project(p, [2], ProjectGenerators()), Ap, bp, BitSet())
    if supportselimination(p, BlockElimination())
        inequality_fulltest(eliminate(p, [1], BlockElimination()), Ap, bp, BitSet())
        inequality_fulltest(project(p, [2], BlockElimination()), Ap, bp, BitSet())
    end
    @test_throws DimensionMismatch project(p, ones(1, 1))
    @test_throws DimensionMismatch project(p, ones(2, 3))
    inequality_fulltest(project(p, reshape([0, 1], 2, 1)), Ap, bp, BitSet())
    inequality_fulltest(fixandeliminate(p, 2, eltype(p) <: AbstractFloat ? 1/4 : 1//4), reshape([1, -1], 2, 1), [1/4, 0], BitSet())
    # tries to create CDDLib.Polyhedron of BigFloat, needs to let CDD decide of the promotion
    #inequality_fulltest(fixandeliminate(p, 2, 1/4, reshape([1, -1], 2, 1), [1/4, 0], BitSet())
    inequality_fulltest(translate(p, -[0, 1]), [1 1; 1 -1;-1 0], [0, 1, 0], BitSet())
    # Testing 2D plotting
    RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)[1].args == ([0, 0.5, 0],
                                                                 [1, 0.5, 0])
end
