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

    @testset "Plots.jl attribute forwarding" begin
        # Creating a simple 2D polyhedron
        v = vrep([0 0; 1 0; 0 1; 1 1])
        p = polyhedron(v)

        # Simulating the user providing custom kwargs to plot()
        # This is exactly what Plots.jl passes into the macro internally
        attrs_with_label = Dict{Symbol, Any}(:label => "Custom Label", :series_annotations => ["1", "2", "3", "4"])
    
        # Calling the recipe directly without needing Plots.jl
        # This returns an array of RecipeData objects
        res_with_label = RecipesBase.apply_recipe(attrs_with_label, p)
        series_attrs = res_with_label[1].plotattributes

        # Verifying custom label prevented `legend --> false` and annotations were kept
        @test !haskey(series_attrs, :legend) || series_attrs[:legend] !== false
        @test haskey(series_attrs, :series_annotations)
        @test series_attrs[:series_annotations] == ["1", "2", "3", "4"]

        # Testing that if the user DOES NOT provide a label, legend defaults to false
        attrs_no_label = Dict{Symbol, Any}()
        res_no_label = RecipesBase.apply_recipe(attrs_no_label, p)
        @test res_no_label[1].plotattributes[:legend] == false
    end
end
