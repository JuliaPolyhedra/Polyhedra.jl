function simplextest(lib::PolyhedraLibrary)
    A = [1 1; -1 0; 0 -1]
    b = [1, 0, 0]
    ls = IntSet([1])
    V = [0 1; 1 0]

    ine = hrep(A, b, ls)
    poly1 = polyhedron(ine, lib)
    @test !isempty(poly1, defaultLPsolverfor(poly1, lpsolver))
    center, radius = chebyshevcenter(poly1, defaultLPsolverfor(poly1, lpsolver))
    @test center ≈ [1/2, 1/2]
    @test radius ≈ 1/2
    inequality_fulltest(poly1, A, b, ls)
    generator_fulltest(poly1, V)

    ext = vrep(V)
    poly2 = polyhedron(ext, lib)
    @test !isempty(poly2)
    inequality_fulltest(poly2, A, b, ls)
    generator_fulltest(poly2, V)

    # x_1 cannot be 2
    poly = polyhedron(hrep([A; 1 0], [b; 2], union(ls, IntSet([4]))), lib)
    @test isempty(poly, defaultLPsolverfor(poly, lpsolver))

    # We now add the vertex (0, 0)
    V0 = [0 0]
    ext0 = vrep(V0)
    @test translate(ext0, [1, 0]).V == [1 0]

    push!(poly1, ext0)
    inequality_fulltest(poly1, A, b, IntSet([]))
    generator_fulltest(poly1, [V0; V])
    push!(poly2, ext0)
    inequality_fulltest(poly2, A, b, IntSet([]))
    generator_fulltest(poly2, [V0; V])

    # nonnegative orthant cut by x_1 + x_2 = 1
    Vray = [1 0; 0 1]
    extray = vrep(zeros(Int, 1, 2), Vray)
    poly3 = polyhedron(extray, lib)
    @test_throws ErrorException chebyshevcenter(poly3)
    Acut = [1 1]
    bcut = [1]
    linsetcut = IntSet([1])
    inecut = hrep(Acut, bcut, linsetcut)
    @test !ininterior([1/2, 1/2], inecut)
    @test inrelativeinterior([1/2, 1/2], inecut)
    push!(poly3, inecut)
    inequality_fulltest(poly3, A, b, ls)
    generator_fulltest(poly3, V)

    # FIXME needs float currently but should be updated
    # poly4 = project(poly1, [1; 0])
    # inequality_fulltest(poly4, [-1; 1], [0, 1], IntSet())
    # generator_fulltest(poly4, [0; 1], [])

    #\
    # \
    # |\
    # |_\
    #    \
    Alin = [1 1]
    blin = [1]
    linsetlin = IntSet(1)
    Vlin = [1 0]
    Rlin = [1 -1]
    inelin = hrep([1 1; -1 -1], [1, -1], IntSet())
    plin = polyhedron(inelin, lib)
    inequality_fulltest(plin, Alin, blin, linsetlin)
    generator_fulltest(plin, Vlin, Rlin, IntSet(), IntSet(1))
    ineout = hrep(plin)
    #@test linset(ineout) == IntSet(1)
    Vlin = [1 0]
    Rlin = [1 -1]
    extlin = vrep(Vlin, [1 -1; -1 1])
    plin = polyhedron(extlin, lib)
    inequality_fulltest(plin, Alin, blin, linsetlin)
    generator_fulltest(plin, Vlin, Rlin, IntSet(), IntSet(1))
    extout = vrep(plin)
    #@test linset(extout) == IntSet(1)
end
