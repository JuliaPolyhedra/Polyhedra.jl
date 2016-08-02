function defaultLPsolverfor{N,T}(p::Rep{N,T})
  if vrepiscomputed(p)
    SimpleVRepPolyhedraModel{N,T}()
  else
    MathProgBase.defaultLPsolver
  end
end
