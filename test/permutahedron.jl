function permutahedrontest{Lib<:PolyhedraLibrary}(lib::Lib)
  A = [1 1 1; 1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]
  b = [6, 3, 3, 3, -1, -1, -1]
  linset = IntSet([1])
  V = [2 3 1; 1 3 2; 3 1 2; 3 2 1; 2 1 3; 1 2 3]
  ine = InequalityDescription(A, b, linset)
  poly = polyhedron(ine, lib)
  @test !isempty(poly)
  inequality_fulltest(poly, A, b, linset)
  generator_fulltest(poly, V, Array(Int, 0, 3))


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
  linsetlift = IntSet([])
  inelift = InequalityDescription(Alift, blift, linsetlift)
  polylift1 = polyhedron(inelift, lib)
  eliminate!(polylift1, IntSet([5,6]))
  polylift0 = eliminate(polylift1)
  inequality_fulltest(polylift0, A, b, linset)
  generator_fulltest(polylift0, V)

# removeredundantinequalities!(polylift0)
# inelift0d = Polyhedra.InequalityDescription{Int}(Polyhedra.Description(inelift0))
# inelift0df = Polyhedra.InequalityDescription{Int}(Polyhedra.Description(inelift0f))
# @test inelift0d.linset == IntSet([1])
# @test length(inelift0d.b) == 7
# @test inelift0d.b[1] / sign(inelift0d.b[1]) == 6
# @test vec(Array{Int}(inelift0d.A[1,:] / sign(inelift0d.b[1]))) == [1; 1; 1] # Array{Int} cast and vec are for julia 0.4
# @test inelift0df.linset == IntSet([1])
# @test length(inelift0df.b) == 7
# @test inelift0df.b[1] / sign(inelift0df.b[1]) == 6
# @test vec(Array{Int}(inelift0df.A[1,:] / sign(inelift0df.b[1]))) == [1; 1; 1] # Array{Int} cast and vec are for julia 0.4
# polylift = CDDPolyhedra(inelift0)
# polyliftf = CDDPolyhedra(inelift0f)
# extunlift = Polyhedra.GeneratorDescription{Int}(Polyhedra.Description(copygenerators(polylift)))
# extunliftf = Polyhedra.GeneratorDescription{Int}(round(Polyhedra.Description(copygenerators(polyliftf))))
# generator_fulltest(extunlift, V, Array(Int, 0, 3))
# generator_fulltest(extunliftf, V, Array(Int, 0, 3))
end
