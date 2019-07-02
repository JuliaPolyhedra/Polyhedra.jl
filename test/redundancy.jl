using Test

using StaticArrays

using Polyhedra

include("solvers.jl")

@testset "Redundancy removal" begin
    @testset "LP-based" begin
        @testset "VRepresentation with $(typeof(x))" for (x, y, z) in [
            ([1, 0], [0, 1], [1, 1]),
            ([1.0, 0.0], [0.0, 1.0], [1.0, 1.0]),
            ((@SVector [1, 0]), (@SVector [0, 1]), (@SVector [1, 1]))]
            T = eltype(x)
            AT = typeof(x)
            d = Polyhedra.FullDim(x)
            @testset "VLinearSpace" begin
                vr = vrep(typeof(Line(x))[]; d=d)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa LinesHull{T, AT}
                @test !haspoints(rm)
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(Line(x))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa LinesHull{T, AT}
                @test npoints(rm) == 1
                @test !hasrays(rm)
                @test nlines(rm) == 1
                vr = convexhull(Line(x), Line(y))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa LinesHull{T, AT}
                @test npoints(rm) == 1
                @test !hasrays(rm)
                @test nlines(rm) == 2
                vr = convexhull(Line(x), Line(y), Line(z))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa LinesHull{T, AT}
                @test npoints(rm) == 1
                @test !hasrays(rm)
                @test nlines(rm) == 2
                vrepm1 = Polyhedra.emptyspace(vr)
                @test vrepm1 === @inferred removevredundancy(vrepm1, lp_solver)
            end
            @testset "Polyhedra.RaysHull" begin
                vr = vrep(typeof(Line(x))[], typeof(Ray(x))[]; d=d)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test !haspoints(rm)
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(Ray(x))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 1
                @test !haslines(rm)
                vr = convexhull(Ray(x), Ray(y))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 2
                @test !haslines(rm)
                for vr in [convexhull(Ray(x), Ray(y), Ray(z)),
                           convexhull(Ray(x), Ray(z), Ray(y)),
                           convexhull(Ray(z), Ray(x), Ray(y))]
                    rm = removevredundancy(vr, lp_solver)
                    @test rm isa Polyhedra.RaysHull{T, AT}
                    @test npoints(rm) == 1
                    @test nrays(rm) == 2
                    rs = collect(rays(rm))
                    @test rs[1] == Ray(x)
                    @test rs[2] == Ray(y)
                    @test !haslines(rm)
                end
                vr = convexhull(Ray(x), Ray(y), Ray(-z))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 3
                @test !haslines(rm)
                vr = convexhull(Ray(x), Ray(y), Line(z))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 2
                @test nlines(rm) == 1
                vr = convexhull(Ray(x), Ray(-y), Line(z))
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 1
                @test nlines(rm) == 1
            end
            @testset "Polyhedra.PointsHull" begin
                vr = vrep(typeof(x)[]; d=d)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.PointsHull{T, AT}
                @test !haspoints(rm)
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(x)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.PointsHull{T, AT}
                @test npoints(rm) == 1
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(x, y)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.PointsHull{T, AT}
                @test npoints(rm) == 2
                @test !hasrays(rm)
                @test !haslines(rm)
                for vr in [convexhull(2x, 2y, z),
                           convexhull(2x, z, 2y),
                           convexhull(z, 2x, 2y)]
                    rm = removevredundancy(vr, lp_solver)
                    @test rm isa Polyhedra.PointsHull{T, AT}
                    @test npoints(rm) == 2
                    ps = collect(points(rm))
                    @test ps[1] == 2x
                    @test ps[2] == 2y
                    @test !hasrays(rm)
                    @test !haslines(rm)
                end
                vr = convexhull(x, y, z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.PointsHull{T, AT}
                @test npoints(rm) == 3
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(x, y, -z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.PointsHull{T, AT}
                @test npoints(rm) == 3
                @test !hasrays(rm)
                @test !haslines(rm)
                vr = convexhull(x, y) + Ray(z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.Hull{T, AT}
                @test npoints(rm) == 2
                @test nrays(rm) == 1
                @test !haslines(rm)
                vr = convexhull(x, y) + Line(z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.Hull{T, AT}
                @test npoints(rm) == 2
                @test !hasrays(rm)
                @test nlines(rm) == 1
                vr = convexhull(x, -y) + Ray(z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.Hull{T, AT}
                @test npoints(rm) == 1
                @test first(points(rm)) == -y
                @test nrays(rm) == 1
                @test !haslines(rm)
                vr = convexhull(x, -y) + Line(z)
                rm = removevredundancy(vr, lp_solver)
                @test rm isa Polyhedra.Hull{T, AT}
                @test npoints(rm) == 1
                @test !hasrays(rm)
                @test nlines(rm) == 1
            end
        end
    end
    @testset "Representation-based" begin
        @testset "isredundant" begin
            @test isredundant(HyperPlane([0, 1], 0) ∩ HalfSpace([1, 1], 1), [0, 0], nl=0)
            @test !isredundant(HalfSpace([1, 1], 1) ∩ HyperPlane([0, 1], 0) ∩ HalfSpace([-1, -1], 1), [1, 0], nl=0)
            #@test !ishredundant(convexhull(SymPoint([1, 0])), HyperPlane([0, 1], 0), d=1)
            #@test isredundant(convexhull(SymPoint([1, 0])), HalfSpace([0, 1], 1), d=1)
            @test isredundant(HalfSpace([1, 1], 2) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0), [1, 1], nl=0)
            @test !isredundant(HalfSpace([1, 1], 2) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0), [1, 1], nl=0, strongly=true)
        end
#    @testset "Non-redundant SymPoint" begin
#        vr = convexhull(SymPoint([-1, 0]), SymPoint([0, 1]))
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1) ∩ HalfSpace([1, 1], 1) ∩ HalfSpace([-1, 1], 1)
#        vr = @inferred Polyhedra.removevredundancy(vr, hr)
#        @test vr isa Polyhedra.PointsHull{2,Int,Vector{Int}}
#        @test collect(sympoints(vr)) == [SymPoint([-1, 0]), SymPoint([0, 1])]
#        @test !haspoints(vr)
#        @test !hasallrays(vr)
#    end
        @testset "Point" begin
            # x - y ≤  1
            # x + y ≥ -1
            hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
            # [-1, 0] is weakly redundant as it belongs to the facet x + y = -1
            @test isredundant(hr, [-1, 0], nl=0)
            @test !isredundant(hr, [-1, 0], strongly=true, nl=0)
            # [0, -1] is not redundant as it is an extreme point
            @test !isredundant(hr, [0, -1], nl=0)
            @test !isredundant(hr, [0, -1], strongly=true, nl=0)
            # [0, 0] is strongly redundant as it is in the relative interior
            @test isredundant(hr, [0, 0], nl=0)
            @test isredundant(hr, [0, 0], strongly=true, nl=0)
            vr = convexhull([-1, 0], [0, -1], [0, 0]) + conichull([1, 1], [-1, 1])
            @test collect(removevredundancy(points(vr), hr, nl=0)) == [[0, -1]]
            @test collect(removevredundancy(points(vr), hr, strongly=true, nl=0)) == [[-1, 0], [0, -1]]
            vrr = Polyhedra.removevredundancy(vr, hr)
            @test collect(points(vrr)) == [[0, -1]]
            @test !haslines(vrr)
            @test collect(rays(vrr)) == [Ray([1, 1]), Ray([-1, 1])]
            vrr = Polyhedra.removevredundancy(vr, hr, strongly=true)
            @test collect(points(vrr)) == [[-1, 0], [0, -1]]
            @test !haslines(vrr)
            @test collect(rays(vrr)) == [Ray([1, 1]), Ray([-1, 1])]
            p = polyhedron(vr)
            Polyhedra.computehrep!(p)
            removevredundancy!(p, strongly=true)
            @test collect(points(p)) == [[-1, 0], [0, -1]]
            @test !haslines(p)
            @test collect(rays(p)) == [Ray([1, 1]), Ray([-1, 1])]
            removevredundancy!(p)
            @test collect(points(p)) == [[0, -1]]
            @test !haslines(p)
            @test collect(rays(p)) == [Ray([1, 1]), Ray([-1, 1])]
        end
#    @testset "Split SymPoint" begin
#        vr = convexhull(SymPoint([-1, 0]), SymPoint([0, 1])) + conichull([1, 1], [-1, 1])
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1)
#        vr = Polyhedra.removevredundancy(vr, hr)
#        @test !hassympoints(vr)
#        @test collect(points(vr)) == [[0, -1]]
#        @test !haslines(vr)
#        @test collect(rays(vr)) == [Ray([1, 1]), Ray([-1, 1])]
#
#        vr = convexhull(SymPoint([1, 0]), SymPoint([0, 1])) + conichull([-1, 1])
#        hr = HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1) ∩ HalfSpace([1, 1], 1)
#        vr = Polyhedra.removevredundancy(vr, hr)
#        @test !hassympoints(vr)
#        @test collect(points(vr)) == [[1, 0], [0, -1]]
#        @test !haslines(vr)
#        @test collect(rays(vr)) == [Ray([-1, 1])]
#    end
    end
end
@testset "Duplicate removal" begin
    h = removeduplicates(hrep([1 1; -1 -1; 1 0; 0 -1], [1, -1, 1, 0]))
    @test nhyperplanes(h) == 1
    @test first(hyperplanes(h)) == HyperPlane([1, 1], 1)
    @test first(hyperplanes(h)) == HyperPlane([-1, -1], -1)
    @test nhalfspaces(h) == 1
    h = removeduplicates(hrep([1 1; -1 -1], [2, -1]))
    @test !hashyperplanes(h)
    @test nhalfspaces(h) == 2
#    v = removeduplicates(convexhull([1, 2], SymPoint([1, 2]), [-1, 1]))
#    @test typeof(v) == Polyhedra.PointsHull{2,Int,Vector{Int}}
#    @test collect(sympoints(v)) == [SymPoint([1, 2])]
#    @test collect(points(v)) == [[-1, 1]]
#    @test !hasallrays(v)
    for v in (#removeduplicates(convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1]) + Line([0, 1])),
              #removeduplicates(Line([0, 1]) + convexhull([1, 2], [2, 1], SymPoint([1, 2]), [-1, 1])),
              removeduplicates(Line([0, 1]) + convexhull([-1, -2], [2, 1], [1, 2], [-1, 1], [-1, -2])),)
        @test typeof(v) == Polyhedra.Hull{Int, Vector{Int}, Int}
        #@test collect(sympoints(v)) == [SymPoint([1, 2])]
        @test collect(points(v)) == [[-1, -2], [2, 1], [1, 2]]
        @test !hasrays(v)
        @test collect(lines(v)) == [Line([0, 1])]
    end
end
