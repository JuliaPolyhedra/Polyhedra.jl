import RecipesBase
function recipetest(lib::Polyhedra.Library)
    v = convexhull([0.0, 0.0], [1.0, 0.0], [0.0, 1.0])
    p = polyhedron(v)
    @test RecipesBase.apply_recipe(Dict{Symbol, Any}(), p)[1].args == ([0.0, 0.0, 1.0, 0.0],
                                                                       [0.0, 1.0, 0.0, 0.0])
    @testset "Error for unbounded polyhedron" begin
        vr = convexhull(v, Ray([1.0, 1.0]))
        pr = polyhedron(vr)
        @test_throws ErrorException RecipesBase.apply_recipe(Dict{Symbol, Any}(), pr)
        vl = convexhull(v, Line([1.0, 1.0]))
        pl = polyhedron(vl)
        @test_throws ErrorException RecipesBase.apply_recipe(Dict{Symbol, Any}(), pl)
    end
end
