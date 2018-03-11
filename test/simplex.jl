function simplextest(lib::PolyhedraLibrary)
    hsim = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)
    vsim = convexhull([0, 1], [1, 0])
    poly1 = polyhedron(hsim, lib)
    @test !isempty(poly1, defaultLPsolverfor(poly1, lpsolver))
    center, radius = chebyshevcenter(poly1, defaultLPsolverfor(poly1, lpsolver))
    @test center ≈ [1/2, 1/2]
    @test radius ≈ 1/2
    @test dim(poly1) == 1 # FIXME doing dim earlier makes chebyshevcenter fail
    inequality_fulltest(poly1, hsim)
    generator_fulltest(poly1, vsim)

    poly2 = polyhedron(vsim, lib)
    @test dim(poly2) == 1
    @test !isempty(poly2)
    inequality_fulltest(poly2, hsim)
    generator_fulltest(poly2, vsim)

    # x_1 cannot be 2
    hempty = hsim ∩ HyperPlane([1, 0], 2)
    poly = polyhedron(hempty, lib)
    @test isempty(poly, defaultLPsolverfor(poly, lpsolver))

    # We now add the vertex (0, 0)
    ext0 = convexhull([0, 0])
    @test collect(points(translate(ext0, [1, 0]))) == [[1, 0]]

    htri = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HalfSpace([1, 1], 1)
    vtri = convexhull(vsim, ext0)
    push!(poly1, ext0)
    inequality_fulltest(poly1, htri)
    generator_fulltest(poly1, vtri)
    push!(poly2, ext0)
    inequality_fulltest(poly2, htri)
    generator_fulltest(poly2, vtri)

    # nonnegative orthant cut by x_1 + x_2 = 1
    extray = conichull(Ray([1, 0]), Ray([0, 1])) # TODO pass the test without "convexhull([0, 0]) + "
    poly3 = polyhedron(extray, lib)
    @test_throws ErrorException chebyshevcenter(poly3)
    @test dim(poly3) == 2
    hcut = intersect(HyperPlane([1, 1], 1))
    vcut = convexhull([1, 0]) + conichull(Line([1, -1]))
    @test !ininterior([1/2, 1/2], hcut)
    @test inrelativeinterior([1/2, 1/2], hcut)
    push!(poly3, hcut)
    @test dim(poly3) == 1
    inequality_fulltest(poly3, hsim)
    generator_fulltest(poly3, vsim)

    # FIXME needs float currently but should be updated
    # poly4 = project(poly1, [1; 0])
    # inequality_fulltest(poly4, [-1; 1], [0, 1], IntSet())
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
    #@test linset(ineout) == IntSet(1)
    vlin = convexhull([1, 0]) + conichull(Ray([1, -1]), Ray([-1, 1]))
    plin = polyhedron(vlin, lib)
    inequality_fulltest(plin, hcut)
    generator_fulltest(plin, vcut)
    #extout = vrep(plin)
    #@test linset(extout) == IntSet(1)
end
