function simplextest{Poly<:Polyhedron}(::Type{Poly}, args...)
  A = [1 1; -1 0; 0 -1]
  b = [1, 0, 0]
  linset = IntSet([1])
  V = [0 1; 1 0]

  ine = InequalityDescription(A, b, linset)
  poly1 = Poly(ine, args...)
  @test !isempty(poly1)
  inequality_fulltest(poly1, A, b, linset)
  generator_fulltest(poly1, V)

  ext = GeneratorDescription(V)
  poly2 = Poly(ext, args...)
  @test !isempty(poly2)
  inequality_fulltest(poly2, A, b, linset)
  generator_fulltest(poly2, V)

  # x_1 cannot be 2
  @test isempty(Poly(Polyhedra.InequalityDescription([A; 1 0], [b; 2], union(linset, IntSet([4]))), args...))

  # We now add the vertex (0, 0)
  V0 = [0 0]
  vertex0 = IntSet([1])
  ext0 = Polyhedra.GeneratorDescription(V0, vertex0)

  push!(poly1, ext0)
  inequality_fulltest(poly1, A, b, IntSet([]))
  generator_fulltest(poly1, [V0; V])
  push!(poly2, ext0)
  inequality_fulltest(poly2, A, b, IntSet([]))
  generator_fulltest(poly2, [V0; V])

  # nonnegative orthant cut by x_1 + x_2 = 1
  Vray = [1 0; 0 1]
  vertexray = IntSet([])
  extray = GeneratorDescription(Vray, vertexray)
  poly3 = Poly(extray, args...)
  Acut = [1 1]
  bcut = [1]
  linsetcut = IntSet([1])
  inecut = Polyhedra.InequalityDescription(Acut, bcut, linsetcut)
  push!(poly3, inecut)
  inequality_fulltest(poly3, A, b, linset)
  generator_fulltest(poly3, V)
end
