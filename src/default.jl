function defaultLPsolverfor{N,T}(p::Rep{N,T})
  if vrepiscomputed(p)
    SimpleVRepSolver()
  else
    MathProgBase.defaultLPsolver
  end
end
