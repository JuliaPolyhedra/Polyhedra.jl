function boardtest{Lib<:PolyhedraLibrary}(lib::Lib)
  A1 = -eye(Int, 9) # x >= 0
  b1 = zeros(Int, 9)
  A2 = eye(Int, 9) # x <= 1
  b2 = ones(Int, 9)
  A3 = zeros(Int, 9, 9)
  b3 = 3 * ones(Int, 9)
  i = 1
  for a = 1:3
    for b = (a+1):3
      for c = 1:3
        for d = (c+1):3
          ac = a + (c-1) * 3
          ad = a + (d-1) * 3
          bc = b + (c-1) * 3
          bd = b + (d-1) * 3
          A3[i, ac] = 1
          A3[i, ad] = 1
          A3[i, bc] = 1
          A3[i, bd] = 1
          i += 1
        end
      end
    end
  end
  A = [A1; A2; A3]
  b = [b1; b2; b3]
  ine = SimpleHRepresentation(A, b)
  poly = polyhedron(ine, lib)
  @fact nrays(poly) --> 0
  @fact npoints(poly) --> 614
  @fact isempty(poly) --> false
  ext  = SimpleVRepresentation(getvrep(poly))
  target = ones(Int, 9) * (3 // 4)
  ok = false
  for i = 1:size(ext.V, 1)
    # In julia v0.4 [i,:] returns a row matrix and in v0.5 it is
    # a 1D vector hence the use of vec
    if vec(ext.V[i,:]) == target
      ok = true
    end
  end
  @fact ok --> true

  cutA = ones(Int, 1, 9)
  cutb = 6
  Acut = [cutA; A]
  bcut = [cutb; b]
  inecut = SimpleHRepresentation(Acut, bcut)
  polycut = polyhedron(inecut, lib)
  @fact isempty(polycut) --> false
  #(isredundant, certificate) = isredundantinequality(polycut, 1)
  #@test !isredundant
  #@test certificate == target
  @fact ishredundant(polycut, 1) --> false
# @test IntSet() == gethredundantindices(polycut) # TODO reactivate it when I figure out why it makes LRS tests fail
 #(issredundant, scertificate) = isstronglyredundantinequality(polycut, 1)
 #@test !issredundant
 #@test scertificate == target
 #@test IntSet([]) == getstronglyredundantinequalities(polycut)
end
