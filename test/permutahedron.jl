function permutahedrontest(lib::Polyhedra.Library)
    A = [1 1 1; 1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]
    b = [6, 3, 3, 3, -1, -1, -1]
    linset = BitSet([1])
    V = [2 3 1; 1 3 2; 3 1 2; 3 2 1; 2 1 3; 1 2 3]
    ine = hrep(A, b, linset)
    poly = polyhedron(ine, lib)
    @test !isempty(poly)
    inequality_fulltest(poly, A, b, linset)
    generator_fulltest(poly, V, Matrix{Int}(undef, 0, 3))

    isa(lib, Polyhedra.Library) && return
    # x1___x4____________1
    #      |         |
    #      V         V
    # x2___x5___x6_______2
    #           |
    #           V
    # x3_________________3
    Alift = [-1  0  0  1  0  0;
              0 -1  0  1  0  0;
             -1 -1  0  1  1  0;
              1  1  0 -1 -1  0;
              0  0 -1  0  0  1;
              0  0  0  0 -1  1;
              0  0 -1  0 -1  1;
              0  0  1  0  1 -1;
              0  0  0 -1  0  0;
              0  0  0  0  0 -1;
              0  0  0 -1  0 -1;
              0  0  0  1  0  1]
    blift = [0; 0; 0; 0; 0; 0; -3; 3; -1; -1; -(1+2); (1+2)]
    linsetlift = BitSet()
    inelift = hrep(Alift, blift, linsetlift)
    polylift2 = polyhedron(inelift, lib)
    #polylift1 = project(polylift2, 1:4)
    polylift1 = eliminate(polylift2, BitSet([5, 6]))
    polylift0 = eliminate(polylift1)
    inequality_fulltest(polylift0, A, b, linset)
    generator_fulltest(polylift0, V)

    # removeredundantinequalities!(polylift0)
    # inelift0d = Polyhedra.HRepresentation{Int}(Polyhedra.Representation(inelift0))
    # inelift0df = Polyhedra.HRepresentation{Int}(Polyhedra.Representation(inelift0f))
    # @test inelift0d.linset == BitSet([1])
    # @test length(inelift0d.b) == 7
    # @test inelift0d.b[1] / sign(inelift0d.b[1]) == 6
    # @test vec(Array{Int}(inelift0d.A[1,:] / sign(inelift0d.b[1]))) == [1; 1; 1] # Array{Int} cast and vec are for julia 0.4
    # @test inelift0df.linset == BitSet([1])
    # @test length(inelift0df.b) == 7
    # @test inelift0df.b[1] / sign(inelift0df.b[1]) == 6
    # @test vec(Array{Int}(inelift0df.A[1,:] / sign(inelift0df.b[1]))) == [1; 1; 1] # Array{Int} cast and vec are for julia 0.4
    # polylift = CDDPolyhedra(inelift0)
    # polyliftf = CDDPolyhedra(inelift0f)
    # extunlift = Polyhedra.VRepresentation{Int}(Polyhedra.Representation(copygenerators(polylift)))
    # extunliftf = Polyhedra.VRepresentation{Int}(round(Polyhedra.Representation(copygenerators(polyliftf))))
    # generator_fulltest(extunlift, V, Array(Int, 0, 3))
    # generator_fulltest(extunliftf, V, Array(Int, 0, 3))
end
