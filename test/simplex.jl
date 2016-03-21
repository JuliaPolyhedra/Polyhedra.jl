function simplextest{Lib<:PolyhedraLibrary}(lib::Lib)
  A = [1 1; -1 0; 0 -1]
  b = [1, 0, 0]
  linset = IntSet([1])
  V = [0 1; 1 0]

  ine = HRepresentation(A, b, linset)
  poly1 = polyhedron(ine, lib)
  @test !isempty(poly1)
  inequality_fulltest(poly1, A, b, linset)
  generator_fulltest(poly1, V)

  ext = VRepresentation(V)
  poly2 = polyhedron(ext, lib)
  @test !isempty(poly2)
  inequality_fulltest(poly2, A, b, linset)
  generator_fulltest(poly2, V)

  # x_1 cannot be 2
  @test isempty(polyhedron(Polyhedra.HRepresentation([A; 1 0], [b; 2], union(linset, IntSet([4]))), lib))

  # We now add the vertex (0, 0)
  V0 = [0 0]
  vertex0 = IntSet([1])
  ext0 = Polyhedra.VRepresentation(V0, vertex0)

  push!(poly1, ext0)
  inequality_fulltest(poly1, A, b, IntSet([]))
  generator_fulltest(poly1, [V0; V])
  push!(poly2, ext0)
  inequality_fulltest(poly2, A, b, IntSet([]))
  generator_fulltest(poly2, [V0; V])

  # nonnegative orthant cut by x_1 + x_2 = 1
  Vray = [1 0; 0 1]
  vertexray = IntSet([])
  extray = VRepresentation(Vray, vertexray)
  poly3 = polyhedron(extray, lib)
  Acut = [1 1]
  bcut = [1]
  linsetcut = IntSet([1])
  inecut = HRepresentation(Acut, bcut, linsetcut)
  push!(poly3, inecut)
  inequality_fulltest(poly3, A, b, linset)
  generator_fulltest(poly3, V)

  poly4 = project(poly1, [1; 0])
  inequality_fulltest(poly4, [-1; 1], [0, 1], IntSet([]))
  generator_fulltest(poly4, [0; 1], [])
end
