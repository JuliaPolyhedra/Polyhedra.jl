import MathProgBase
function simplextest(lib::Polyhedra.Library)
    hsim = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)
    vsim = convexhull([0, 1], [1, 0])
    poly1 = polyhedron(hsim, lib)
    @test !isempty(poly1)
    center, radius = chebyshevcenter(poly1)
    @test center ≈ [1/2, 1/2]
    @test radius ≈ 1/2
    @test dim(poly1) == 1 # FIXME doing dim earlier makes chebyshevcenter fail
    inequality_fulltest(poly1, hsim)
    generator_fulltest(poly1, vsim)

    # Test incidence
    hpidx = first(eachindex(hyperplanes(poly1)))
    for ps in (incidentpoints(poly1, hpidx), get.(poly1, incidentpointindices(poly1, hpidx)))
        @test (ps == [[0, 1], [1, 0]] || ps == [[1, 0], [0, 1]])
    end
    for hsidx in eachindex(halfspaces(poly1))
        h = get(poly1, hsidx)
        if dot(h.a, [1, 0]) ≈ h.β
            expps = [[1, 0]]
        else
            @assert dot(h.a, [0, 1]) ≈ h.β
            expps = [[0, 1]]
        end
        for ps in (incidentpoints(poly1, hsidx), get.(poly1, incidentpointindices(poly1, hsidx)))
            @test ps == expps
        end
    end
    for pidx in eachindex(points(poly1))
        for hps in (incidenthyperplanes(poly1, pidx), get.(poly1, incidenthyperplaneindices(poly1, pidx)))
            @test hps == [HyperPlane([1, 1], 1)]
        end
        for hss in (incidenthalfspaces(poly1, pidx), get.(poly1, incidenthalfspaceindices(poly1, pidx)))
            h = hss[1]
            @test dot(h.a, get(poly1, pidx)) ≈ h.β
        end
    end

    @test_throws DimensionMismatch MathProgBase.linprog(ones(3), poly1)
    sol = MathProgBase.linprog([-2, 0], poly1)
    @test sol.status == :Optimal
    @test sol.objval == -2
    @test sol.sol == [1, 0]
    sol = MathProgBase.linprog([-1, -3], poly1)
    @test sol.status == :Optimal
    @test sol.objval == -3
    @test sol.sol == [0, 1]

    poly2 = polyhedron(vsim, lib)
    @test dim(poly2) == 1
    @test !isempty(poly2)
    inequality_fulltest(poly2, hsim)
    generator_fulltest(poly2, vsim)

    # x_1 cannot be 2
    hempty = hsim ∩ HyperPlane([1, 0], 2)
    poly = polyhedron(hempty, lib)
    @test isempty(poly)

    # We now add the vertex (0, 0)
    ext0 = convexhull([0, 0])
    @test collect(points(translate(ext0, [1, 0]))) == [[1, 0]]

    htri = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HalfSpace([1, 1], 1)
    vtri = convexhull(vsim, ext0)
    convexhull!(poly1, ext0)
    inequality_fulltest(poly1, htri)
    generator_fulltest(poly1, vtri)
    convexhull!(poly2, ext0)
    inequality_fulltest(poly2, htri)
    generator_fulltest(poly2, vtri)

    # nonnegative orthant cut by x_1 + x_2 = 1
    vray = conichull([1, 0], [0, 1])
    poly3 = polyhedron(vray, lib)
    @test_throws ErrorException chebyshevcenter(poly3)
    @test dim(poly3) == 2

    @testset "LinProg" begin
        sol = MathProgBase.linprog([1, 1], poly3)
        @test sol.status == :Optimal
        @test sol.objval == 0
        @test sol.sol == [0, 0]
        sol = MathProgBase.linprog([0, -1], poly3)
        @test sol.status == :Unbounded
        @test sol.attrs[:unboundedray] == [0, 1]
    end

    hcutel = HyperPlane([1, 1], 1)
    hcut = intersect(hcutel)
    vcut = convexhull([1, 0]) + conichull(Line([1, -1]))
    @test !ininterior([1/2, 1/2], hcut)
    @test inrelativeinterior([1/2, 1/2], hcut)

    poly4 = copy(poly3)

    polycut3 = poly3 ∩ hcutel
    @test dim(polycut3) == 1
    inequality_fulltest(polycut3, hsim)
    generator_fulltest(polycut3, vsim)
    intersect!(poly3, hcutel)
    @test dim(poly3) == 1
    inequality_fulltest(poly3, hsim)
    generator_fulltest(poly3, vsim)

    # It should not have been cut as it is a copy of poly3
    @test dim(poly4) == 2

    polycut4 = poly4 ∩ hcut
    @test dim(polycut4) == 1
    inequality_fulltest(polycut4, hsim)
    generator_fulltest(polycut4, vsim)
    intersect!(poly4, hcut)
    @test dim(poly4) == 1
    inequality_fulltest(poly4, hsim)
    generator_fulltest(poly4, vsim)

    # FIXME needs float currently but should be updated
    # poly4 = project(poly1, [1; 0])
    # inequality_fulltest(poly4, [-1; 1], [0, 1], BitSet())
    # generator_fulltest(poly4, [0; 1], [])

    #\
    # \
    # |\
    # |_\
    #    \
    hlin = HalfSpace([1, 1], 1) ∩ HalfSpace([-1, -1], -1)
    plin = polyhedron(hlin, lib)
    @test dim(plin) == 1
    inequality_fulltest(plin, hcut)
    generator_fulltest(plin, vcut)
    #ineout = hrep(plin)
    #@test linset(ineout) == BitSet(1)
    vlin = convexhull([1, 0]) + conichull([1, -1], [-1, 1])
    plin = polyhedron(vlin, lib)
    inequality_fulltest(plin, hcut)
    generator_fulltest(plin, vcut)
    #extout = vrep(plin)
    #@test linset(extout) == BitSet(1)
end
