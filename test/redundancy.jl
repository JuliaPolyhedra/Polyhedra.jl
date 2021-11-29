using Test

using StaticArrays

using Polyhedra

include("solvers.jl")
include("dummy_optimizer.jl")

@testset "Redundancy removal" begin
    @testset "LP-based" begin
        @testset "HRepresentation with $(typeof(x))" for (x, y, z, w) in [
            ([1, 0], [0, 1], [1, 1], [0, 0]),
            ([1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]),
            ((@SVector [1, 0]), (@SVector [0, 1]), (@SVector [1, 1]), (@SVector [0, 0]))]
            T = eltype(x)
            AT = typeof(x)
            d = Polyhedra.FullDim(x)
            @testset "HAffineSpace" begin
                for hr in [hrep(typeof(HyperPlane(x, one(T)))[]; d=d),
                           intersect(HyperPlane(w, zero(T)))]
                    @test dim(hr) == 2
                    for h in [detecthlinearity(hr, lp_solver),
                              removehredundancy(hr, lp_solver)]
                        @test dim(h) == 2
                        @test h isa HyperPlanesIntersection{T, AT}
                        @test !hashyperplanes(h)
                        @test !hashalfspaces(h)
                    end
                end
                hr = intersect(HyperPlane(x, one(T)))
                @test dim(hr) == 1
                for h in [detecthlinearity(hr, lp_solver),
                          removehredundancy(hr, lp_solver)]
                    @test dim(h) == 1
                    @test h isa HyperPlanesIntersection{T, AT}
                    @test nhyperplanes(h) == 1
                    @test !hashalfspaces(h)
                end
                hr = intersect(HyperPlane(w, one(T)))
                @test dim(hr) == -1
                for h in [detecthlinearity(hr, lp_solver),
                          removehredundancy(hr, lp_solver)]
                    @test dim(h) == -1
                    @test h isa HyperPlanesIntersection{T, AT}
                    @test nhyperplanes(h) == 3
                    @test !hashalfspaces(h)
                end
            end
            @testset "Intersection" begin
                for hr in [HalfSpace(z + x, 3) ∩ HalfSpace(-x, -one(T)) ∩ HalfSpace(-2y, -2),
                           HalfSpace(x, one(T)) ∩ HalfSpace(y, one(T)) ∩ HalfSpace(-z, -2)]
                    h = detecthlinearity(hr, lp_solver, verbose=2)
                    @test h isa Polyhedra.Intersection{T, AT}
                    @test nhyperplanes(h) == 2
                    @test !hashalfspaces(h)
                end
                for hr in [HalfSpace(z, zero(T)) ∩ HalfSpace(-z, -zero(T)),
                           HalfSpace(x, one(T)) ∩ HalfSpace(-x, -one(T))]
                    h = detecthlinearity(hr, lp_solver, verbose=2)
                    @test h isa Polyhedra.Intersection{T, AT}
                    @test nhyperplanes(h) == 1
                    @test !hashalfspaces(h)
                end
                for hr in [HalfSpace(-x, one(T)) ∩ HalfSpace(-y, one(T)) ∩ HalfSpace(x, one(T)) ∩ HalfSpace(y, one(T)) ∩ HalfSpace(-x, -9),
                           HalfSpace(x, zero(T)) ∩ HalfSpace(-x, -one(T))]
                    h = detecthlinearity(hr, lp_solver, verbose=2)
                    @test h isa Polyhedra.Intersection{T, AT}
                    @test nhyperplanes(h) == 3
                    @test !hashalfspaces(h)
                end
            end
        end
        @testset "VRepresentation with $(typeof(x))" for (x, y, z, w) in [
            ([1, 0], [0, 1], [1, 1], [0, 0]),
            ([1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]),
            ((@SVector [1, 0]), (@SVector [0, 1]), (@SVector [1, 1]), (@SVector [0, 0]))]
            T = eltype(x)
            AT = typeof(x)
            d = Polyhedra.FullDim(x)
            @testset "VLinearSpace" begin
                for vr in [vrep(typeof(Line(x))[]; d=d),
                           convexhull(Line(w))]
                    for v in [detectvlinearity(vr, lp_solver),
                              removevredundancy(vr, lp_solver)]
                        @test v isa LinesHull{T, AT}
                        @test !haspoints(v)
                        @test !hasrays(v)
                        @test !haslines(v)
                    end
                end
                vr = convexhull(Line(x))
                for v in [detectvlinearity(vr, lp_solver),
                          removevredundancy(vr, lp_solver)]
                    @test v isa LinesHull{T, AT}
                    @test npoints(v) == 1
                    @test !hasrays(v)
                    @test nlines(v) == 1
                end
                for vr in [convexhull(Line(x), Line(y)),
                           convexhull(Line(x), Line(y), Line(z)),
                           convexhull(Line(x), Line(z), Line(y)),
                           convexhull(Line(x), Line(z), Line(z))]
                    for v in [detectvlinearity(vr, lp_solver),
                              removevredundancy(vr, lp_solver)]
                        @test v isa LinesHull{T, AT}
                        @test npoints(v) == 1
                        @test !hasrays(v)
                        @test nlines(v) == 2
                    end
                end
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
                vr = conichull(Ray(x), Ray(y), Ray(-z))
                rm = detectvlinearity(vr, lp_solver)
                @test rm isa Polyhedra.RaysHull{T, AT}
                @test npoints(rm) == 1
                @test nrays(rm) == 0
                @test nlines(rm) == 2
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
            p = polyhedron(vr, DefaultLibrary{Rational{BigInt}}(Polyhedra.OppositeMockOptimizer))
            Polyhedra.computehrep!(p)
            removevredundancy!(p, strongly=true)
            @test collect(points(p)) == [[-1, 0], [0, -1]]
            @test !haslines(p)
            @test collect(rays(p)) == [Ray([1, 1]), Ray([-1, 1])]
            removevredundancy!(p, planar=false)
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
        @testset "Solver error in redundancy" begin
            vr = convexhull([1.0, 0.0], [0.0, 1.0])
            err = ErrorException("Solver returned `OTHER_ERROR` when attempting to determine whether `[1.0, 0.0]` of index `Polyhedra.Index{Float64, Vector{Float64}}(1)` is redundant. Solver specific status: Dummy status")
            @test_throws err removevredundancy(vr, DummyOptimizer)
        end
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

@testset "Numerically delicate #270" begin
    A = [
        -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        -0.00106131548486069 0.016634209829441304 -0.0006666769554416961 -1.0 -0.0 -0.0 -0.0 -0.0
        -0.0004116155738065077 0.01693983068952557 -2.6779479610930705e-5 -0.0 -1.0 -0.0 -0.0 -0.0
        -0.0024283354284051484 0.08467041915859919 -0.0005000035065344641 -0.0 -0.0 -1.0 -0.0 -0.0
        -0.002024394302735674 0.08367543404800237 -0.00013321128253368072 -0.0 -0.0 -0.0 -1.0 -0.0
        -0.0009309923368121566 0.009059293773203692 -0.0006675849264721348 -0.0 -0.0 -0.0 -0.0 -1.0
        0.00106131548486069 -0.016634209829441304 0.0006666769554416961 1.0 0.0 0.0 0.0 0.0
        0.0004116155738065077 -0.01693983068952557 2.6779479610930705e-5 0.0 1.0 0.0 0.0 0.0
        0.0024283354284051484 -0.08467041915859919 0.0005000035065344641 0.0 0.0 1.0 0.0 0.0
        0.002024394302735674 -0.08367543404800237 0.00013321128253368072 0.0 0.0 0.0 1.0 0.0
        0.0009309923368121566 -0.009059293773203692 0.0006675849264721348 0.0 0.0 0.0 0.0 1.0
        -0.08320472072273968 2.548630897981053 0.001134866706464998 0.0 0.0 0.0 0.0 0.0
        1210.55555992692 -37126.858723322744 -13.95377106842899 -0.0 -0.0 -0.0 -0.0 -0.0
        -1860.676468633397 57155.38024488611 18.21165578673308 -0.0 -0.0 -0.0 -0.0 -0.0
    ]
    # The last row is almost considered an equality, depends on `ztol`.
    # If it's detected as an equality then the polyhedron is detected as empty.
    # See https://github.com/JuliaPolyhedra/Polyhedra.jl/issues/270
    b = [-20.995, -1.007487377909021, -0.9877340129555137, -0.9705208833525414, -0.9548402631790559, -0.025711479466835108, 1.007487377909021, 0.9877340129555137, 0.9705208833525414, 0.9548402631790559, 0.025711479466835108, 0.8046759919340509, -11751.33121007926, 18148.933864538412]
    lib = DefaultLibrary{Float64}(lp_solver)
    P = polyhedron(hrep(A, b), lib)
    removevredundancy!(P, ztol=1e-4)
    @test nhyperplanes(P) == 7
    @test nhalfspaces(P) == 2
    @test npoints(P) > 0
end

@testset "Noise linearity #259" begin
    h = HalfSpace([-1.3333333333333333, -2.220446049250313e-16], 0.0) ∩ HalfSpace([0.0, 1.333333333333333], 0.0) ∩ HalfSpace([-2.7755575615628914e-17, -0.6666666666666667], 0.0) ∩ HalfSpace([-0.33333333333333337, 0.0], 0.0)
    # Presolve off on purpose to test that `λ_l∞` fixes it.
    hh = detecthlinearity(h, GLPK.Optimizer)
    @test nhyperplanes(hh) == 1
    hhh = removehredundancy(hh, GLPK.Optimizer)
    @test nhyperplanes(hhh) == 1
    @test nhalfspaces(hhh) == 1
end
