using Test
using Polyhedra
using StaticArrays
using LinearAlgebra

function empty_test(h::HRepresentation{T}) where {T}
    #v = @inferred doubledescription(h)
    v = doubledescription(h)
    U = T == Float64 ? Float64 : Rational{BigInt}
    v = doubledescription(h)
    @test v isa Polyhedra.Hull{U,Vector{U}}
    @test !haspoints(v)
    @test !hasallrays(v)
end
empty_hyperplane_test(z) = empty_test(intersect(HyperPlane([z, z], 1)))
empty_halfspace_test(z) = empty_test(intersect(HalfSpace([z, z], -1)))
empty_space_test(z) = empty_test(HyperPlane([z, 1], 1) ∩ HyperPlane([z, 1], -1))
empty_range_test(z) = empty_test(HalfSpace([1, z], -1) ∩ HalfSpace([-1, z], -1))

function simplex_test(z)
    h = hrep([
        HalfSpace([1, 1], 1),
        HalfSpace([1, -1], z),
        HalfSpace([-1, z], z),
    ])
    v = doubledescription(h)
    T = z isa AbstractFloat ? Float64 : Rational{BigInt}
    @test v isa Polyhedra.Hull{T,Vector{T}}
    @test collect(points(v)) == [[z, z], [z, 1], [1//2, 1//2]]
    @test !hasallrays(v)
    h = hrep([HalfSpace(zeros(typeof(z), 2), z); h.halfspaces])
    h = hrep([
        HalfSpace([1, 1], 1),
        HalfSpace([z, z], z),
        HalfSpace([1, -1], z),
        HalfSpace([-1, z], z),
    ])
    #intersect!(h, )
    v = doubledescription(h)
    @test v isa Polyhedra.Hull{T,Vector{T}}
    @test collect(points(v)) == [[z, z], [z, 1], [1//2, 1//2]]
    @test !hasallrays(v)
end

@testset "Double Description" begin
    @test Polyhedra.polytypefor(Float32) == Float64
    @testset "H-representation -> V-representation" begin
        @testset "Intersection" begin
            @testset "Vector" begin
                @testset "Exact" begin
                    simplex_test(0)
                end
                @testset "Numerical" begin
                    simplex_test(0.0)
                end
            end
            @testset "SVector" begin
                @testset "Exact" begin
                    h = hrep([
                        HalfSpace((@SVector [1, 1]), 1),
                        HalfSpace((@SVector [1, -1]), 0),
                        HalfSpace((@SVector [-1, 0]), 0),
                    ])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{
                        Rational{BigInt},
                        SVector{2,Rational{BigInt}},
                    }
                    @test collect(points(v)) == [
                        (@SVector [0, 0]),
                        (@SVector [0, 1]),
                        (@SVector [1//2, 1//2]),
                    ]
                    @test !hasallrays(v)
                end
                @testset "Numerical" begin
                    h = hrep([
                        HalfSpace((@SVector [1.0, 1]), 1),
                        HalfSpace((@SVector [1.0, -1]), 0),
                        HalfSpace((@SVector [-1.0, 0]), 0),
                    ])
                    #v = @inferred doubledescription(h)
                    v = doubledescription(h)
                    @test v isa Polyhedra.Hull{Float64,SVector{2,Float64}}
                    @test collect(points(v)) == [
                        (@SVector [0.0, 0.0]),
                        (@SVector [0.0, 1.0]),
                        (@SVector [1/2, 1/2]),
                    ]
                    @test !hasallrays(v)
                end
            end
        end
        @testset "MixedMatHRep" begin
            @testset "Exact" begin
                h = hrep(
                    [
                        1 1
                        1 -1
                        -1 0
                    ],
                    [1, 0, 0],
                )
                #v = @inferred doubledescription(h)
                v = doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{Rational{BigInt}}
                @test v.V == [0//1 0//1; 0//1 1//1; 1//2 1//2]
                @test v.R == zeros(Rational{BigInt}, 0, 2)
                @test isempty(v.Rlinset)
            end
            @testset "Numerical" begin
                h = hrep(
                    [
                        1 1
                        1 -1
                        -1 0
                    ],
                    [1.0, 0, 0],
                )
                #v = @inferred doubledescription(h)
                v = doubledescription(h)
                @test v isa Polyhedra.MixedMatVRep{Float64}
                @test v.V == [0 0; 0 1; 1/2 1/2]
                @test v.R == zeros(0, 2)
                @test isempty(v.Rlinset)
            end
        end
        _str(::Float64) = "Numerical"
        _str(::Int) = "Exact"
        @testset "Empty Intersection $(_str(z))" for z in [0, 0.0]
            @testset "0x_1 + 0x_2 = 1" begin
                empty_hyperplane_test(z)
            end
            @testset "0x_1 + 0x_2 ≤ -1" begin
                empty_halfspace_test(z)
            end
            @testset "-1 = 0x_1 + x_2 = 1" begin
                empty_space_test(z)
            end
            @testset "1 ≤ x_1 + 0x_2 ≤ -1" begin
                empty_range_test(z)
            end
        end
        @testset "Simple hyperplane $(_str(o))" for o in [1, 1.0]
            h = intersect(HyperPlane([o, o], o))
            v = doubledescription(h)
            @test nlines(v) == 1
            @test first(lines(v)) == Line([-1, 1])
            @test !hasrays(v) == 1
            @test npoints(v) == 1
            @test first(points(v)) == [1, 0]
        end
        @testset "Quadrilateral $(_str(z))" for z in [0, 0.0]
            h =
                HalfSpace([2, 1], 4) ∩ HalfSpace([1, 2], 4) ∩
                HalfSpace([-1, z], z) ∩ HalfSpace([z, -1], z)
            v = doubledescription(h)
            @test !haslines(v)
            @test !hasrays(v)
            @test npoints(v) == 4
            @test all(points(v) .≈ [[0, 0], [2, 0], [0, 2], [4//3, 4//3]])
        end
    end
    @testset "V-representation -> H-representation" begin
        @testset "Dehomogenization ($T)" for T in (Float64, Rational{BigInt})
            v = vrep(zeros(T, 1, 2), T[-1 0; 0 -1])
            h = doubledescription(v)
            @test !hashyperplanes(h)
            @test collect(halfspaces(h)) ==
                  [HalfSpace(T[1, 0], zero(T)), HalfSpace(T[0, 1], zero(T))]
            # Translating the cone preserves genuine nonzero offsets.
            h = doubledescription(vrep(T[2 3], T[-1 0; 0 -1]))
            @test collect(halfspaces(h)) ==
                  [HalfSpace(T[1, 0], T(2)), HalfSpace(T[0, 1], T(3))]
            # Retain equalities for a cone in a proper subspace.
            h = doubledescription(vrep(zeros(T, 1, 2), T[-1 0]))
            @test nhyperplanes(h) == 1
            @test first(hyperplanes(h)).β == 0
            @test abs.(first(hyperplanes(h)).a) == T[0, 1]
            @test collect(halfspaces(h)) == [HalfSpace(T[1, 0], zero(T))]
            # Full space has no inequalities, but still has dimension two.
            h = doubledescription(
                vrep([Line(T[1, 0]), Line(T[0, 1])], Ray{T,Vector{T}}[]),
            )
            @test fulldim(h) == 2
            @test !hashyperplanes(h)
            @test !hashalfspaces(h)
        end
        @testset "Hull" begin
            @testset "Vector" begin
                @testset "Exact" begin
                    v = conichull([1, 0], [0, 1])
                    h = doubledescription(v)
                    #h = @inferred doubledescription(v)
                    @test !hashyperplanes(h)
                    @test collect(halfspaces(h)) ==
                          [HalfSpace([-1, 0], 0), HalfSpace([0, -1], 0)]
                end
                @testset "Numerical" begin
                    v = conichull([1.0, 0.0], [0.0, 1.0])
                    h = doubledescription(v)
                    #h = @inferred doubledescription(v)
                    @test !hashyperplanes(h)
                    @test collect(halfspaces(h)) ==
                          [HalfSpace([-1, 0], 0), HalfSpace([0, -1], 0)]
                end
            end
        end
    end
end

# See https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/251
@testset "Ill-conditioned cut line #251" begin
    V1 = [
        [
            0.26944343620620703,
            0.5820831768940575,
            0.043944752271130684,
            0.0,
            0.0,
            0.0,
            0.0,
            0.002290710865434633,
        ],
        [
            0.269764654943868,
            0.5826024342879536,
            0.04424939599453754,
            -0.3243733749380128,
            -0.27971766698656575,
            -0.040712222052006995,
            -0.4230541656638591,
            0.0020389818225279343,
        ],
        [
            0.269764654943868,
            0.5826024342879536,
            0.04424939599453754,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0020389818225279343,
        ],
        [
            0.3608261278169685,
            0.6795241675058606,
            0.13674690705026338,
            -0.4026336283204254,
            -0.3521265118897102,
            -0.0985458649168468,
            -0.4867600723992247,
            0.022141261211450644,
        ],
        [
            0.3608261278169685,
            0.6795241675058606,
            0.13674690705026338,
            0.0,
            0.0,
            0.0,
            0.0,
            0.022141261211450644,
        ],
        [
            0.37879544708728713,
            0.7586303128811639,
            0.11059986791879892,
            -0.4390178801299712,
            -0.37783324554721326,
            -0.08204843892446032,
            -0.5496916242909842,
            0.017745282399702605,
        ],
        [
            0.37879544708728713,
            0.7586303128811639,
            0.11059986791879892,
            0.0,
            0.0,
            0.0,
            0.0,
            0.017745282399702605,
        ],
        [
            0.35699186847910663,
            0.7349374913159256,
            0.08851175424220448,
            -0.42014699571976977,
            -0.36009900240651904,
            -0.06804164876801723,
            -0.534431302171298,
            0.013816230720892301,
        ],
        [
            0.35699186847910663,
            0.7349374913159256,
            0.08851175424220448,
            0.0,
            0.0,
            0.0,
            0.0,
            0.013816230720892301,
        ],
        [
            0.3566093175846898,
            0.734546246932556,
            0.08812122743510388,
            -0.41982256106796995,
            -0.3598078125821968,
            -0.06780390511187857,
            -0.5341638922672701,
            0.013702811188085457,
        ],
    ]
    V2 = [
        zeros(8),
        [
            0.2697646549,
            0.582602434,
            0.044249395,
            0.3243733749,
            0.279717666,
            0.040712222,
            0.423054165,
            0.00203898,
        ],
        [0.2697646549, 0.582602434, 0.044249395, 0, 0, 0, 0, 0.00203898],
        [
            0.3608261278,
            0.679524167,
            0.136746907,
            0.402633628,
            0.35212651,
            0.098545864,
            0.486760072,
            0.02214126,
        ],
        [0.3608261278, 0.679524167, 0.136746907, 0, 0, 0, 0, 0.02214126],
        [
            0.378795447,
            0.7586303128,
            0.110599867,
            0.43901788,
            0.377833245,
            0.0820484389,
            0.549691624,
            0.017745282,
        ],
        [0.378795447, 0.7586303128, 0.110599867, 0, 0, 0, 0, 0.017745282],
        [
            0.35699186847,
            0.7349374913,
            0.08851175424,
            0.42014699571,
            0.3600990024,
            0.06804164876,
            0.53443130217,
            0.0138162307,
        ],
        [0.35699186847, 0.7349374913, 0.08851175424, 0, 0, 0, 0, 0.0138162307],
        [
            0.35660931758,
            0.7345462469,
            0.08812122743,
            0.419822561067,
            0.35980781258,
            0.06780390511,
            0.534163892267,
            0.0137028111,
        ],
    ]
    for V in [V1, V2]
        h = doubledescription(vrep(V))
        @test !hashyperplanes(h)
        # The exact polytope has 24 facets but some of them are numerically
        # redundant, we only check that the result is consistent with `V`.
        @test all(halfspaces(h)) do hs
            all(v -> dot(hs.a, v) <= hs.β + 1e-6 * maximum(abs, hs.a), V)
        end
    end
end
