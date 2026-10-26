using Test
using Polyhedra
using StaticArrays
using SparseArrays
import Aqua
import GeometryBasics

@testset "Aqua regressions" begin
    @testset "Typed full dimension" begin
        @test Polyhedra.FullDim === Union{Int,StaticArrays.Size}
        @test (@inferred Polyhedra.typed_fulldim(Vector{Int})) == -1
        @test (@inferred Polyhedra.typed_fulldim(SVector{2,Int})) == Size(2)
        for x in ([1, 0], sparsevec([1], [1], 2), SVector(1, 0))
            d = x isa SVector ? Size(2) : 2
            @test (@inferred Polyhedra.typed_fulldim(x)) == d
            for el in (HyperPlane(x, 1), HalfSpace(x, 1), Line(x), Ray(x))
                @test (@inferred Polyhedra.typed_fulldim(el)) == d
                @test fulldim(el) == 2
            end
            for rep in (
                hrep([HyperPlane(x, 1)]),
                hrep([HalfSpace(x, 1)]),
                vrep([x]),
                vrep([Line(x)]),
                vrep([Ray(x)]),
            )
                @test (@inferred Polyhedra.typed_fulldim(rep)) == d
                if !(rep isa Polyhedra.HAffineSpace)
                    @test Polyhedra.typed_fulldim(polyhedron(rep)) == d
                end
                @test fulldim(rep) == 2
            end
        end
    end
    @testset "Representation iterators" begin
        for T in (Int, Float64, Rational{BigInt})
            x, y = T[1, 0], T[0, 1]
            v = vrep([x], [Line(y)], [Ray(x)])
            h = hrep([HyperPlane(x, one(T))], [HalfSpace(y, one(T))])
            for (f, rep) in (
                (points, v),
                (lines, v),
                (rays, v),
                (hyperplanes, h),
                (halfspaces, h),
            )
                expected = collect(f(rep))
                @test collect(@inferred f(rep, rep)) == [expected; expected]
                @test_throws MethodError f()
            end
            @test length(Polyhedra.hreps(h)) == 2
            @test length(Polyhedra.hreps(h, h)) == 2
            affine = hrep([HyperPlane(x, one(T))])
            @test length(Polyhedra.hreps(affine, affine)) == 1
            @test collect(only(Polyhedra.hreps(affine, affine))) ==
                  [HyperPlane(x, one(T)), HyperPlane(x, one(T))]
        end
        # Combining coefficient types requires an explicit conversion, just
        # as for SingleRepIterator and RepIterator.
        @test_throws MethodError points(vrep([[1, 0]]), vrep([[1.0, 0.0]]))
        @test Polyhedra.rreps() == ()
        @test Polyhedra.rmap(identity, 2, Float64) == ()
    end
    @testset "Cone sums and empty representations" begin
        for T in (Int, Float64, Rational{BigInt})
            x, y = T[1, 0], T[0, 1]
            cone = conichull(x)
            other = conichull(y)
            @test Polyhedra.sumpoints(2, T, cone, other) == ()
            total = cone + other
            @test collect(rays(total)) == [Ray(x), Ray(y)]
            @test collect(points(total)) == [zeros(T, 2)]
            space = convexhull(Line(x)) + convexhull(Line(y))
            @test collect(lines(space)) == [Line(x), Line(y)]
            empty = Polyhedra.emptyspace(vrep([x]))
            @test (@inferred removevredundancy(
                empty,
                hrep([HalfSpace(x, one(T))]),
            )) === empty
            @test removevredundancy(
                empty,
                hrep([HalfSpace(x, one(T))]);
                strongly = true,
            ) === empty
        end
    end
    @testset "Mesh initialization" begin
        ext = Base.get_extension(Polyhedra, :PolyhedraGeometryBasicsExt)
        Aqua.test_unbound_args(ext)
        for N in (2, 3), static in (false, true)
            # Build the simplex from the origin and coordinate vectors.
            vertices = hcat(zeros(N), [Float64(i == j) for i in 1:N, j in 1:N])
            ps = [
                static ? SVector{N}(vertices[:, j]) : vertices[:, j] for
                j in axes(vertices, 2)
            ]
            p = polyhedron(vrep(ps))
            mesh = Polyhedra.Mesh(p)
            @test mesh isa ext.Mesh{N,Float64}
            @test isempty(mesh.coordinates)
            @test isempty(mesh.faces)
            @test isempty(mesh.normals)
            @test eltype(mesh.coordinates) == GeometryBasics.Point{N,Float64}
            @test eltype(mesh.faces) == GeometryBasics.TriangleFace{Int}
            @test eltype(mesh.normals) == GeometryBasics.Point{3,Float64}
            coords = GeometryBasics.coordinates(mesh)
            @test !isempty(coords)
            @test GeometryBasics.coordinates(mesh) === coords
            @test length(coords) == 3length(GeometryBasics.faces(mesh))
            @test length(GeometryBasics.normals(mesh)) == length(coords)
        end
    end
end
