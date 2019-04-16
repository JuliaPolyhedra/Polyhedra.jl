import RecipesBase
function recipetest(lib::Polyhedra.Library)
    @testset "Error for 3D polyhedron" begin
        err = ErrorException("Plotting 3-dimensional polyhedron with Plots is not supported, use Makie or MeshCat.")
        v = convexhull(ones(3), zeros(3))
        p = polyhedron(v, lib)
        @test_throws err RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)
    end
    @testset "Error for 1D polyhedron" begin
        err = ErrorException("Plotting 1-dimensional polyhedron with Plots is not supported.")
        h = HalfSpace([1], 1) ∩ HalfSpace([-1], 1)
        p = polyhedron(h, lib)
        @test_throws err RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)
    end
    v = convexhull([0.0, 0.0], [1.0, 0.0], [0.0, 1.0])
    p = polyhedron(v)
    @test RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)[1].args == ([0.0, 0.0, 1.0, 0.0],
                                                                       [0.0, 1.0, 0.0, 0.0])
    @testset "Error for unbounded polyhedron" begin
        err = ErrorException("Rays not supported yet in the 2D plotting recipe.")
        vr = convexhull(v, Ray([1.0, 1.0]))
        pr = polyhedron(vr, lib)
        @test_throws err RecipesBase.apply_recipe(Dict{Symbol, Any}(), pr)
        vl = convexhull(v, Line([1.0, 1.0]))
        pl = polyhedron(vl, lib)
        @test_throws err RecipesBase.apply_recipe(Dict{Symbol, Any}(), pl)
    end
    @testset "Error for empty polyhedron" begin
        err = ErrorException("Plotting empty polyhedron is not supported.")
        h = HalfSpace([1, 1], 0) ∩ HyperPlane([1, 1], 1)
        p = polyhedron(h, lib)
        @test_throws err RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)
    end
end
