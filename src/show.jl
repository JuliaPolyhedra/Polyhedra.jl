import Base.show

function Base.show{N,T}(io::IO, rep::Representation{N,T})
  if typeof(rep) <: HRepresentation
    print(io, "H")
  else
    print(io, "V")
  end
  println(io, "-representation")

  ls = linset(rep)
  if !isempty(ls)
    print(io, "linearity $(length(ls))");
    for i in ls
      print(io, " $i")
    end
    println(io)
  end

  println(io, "begin")
  if T <: AbstractFloat
    typename = "real"
  elseif T <: Integer
    typename = "integer"
  else
    typename = "rational"
  end
  println(io, " $(length(rep)) $(N+1) $typename")
  if typeof(rep) <: HRepresentation
    for h in hrep(rep)
      print(io, " $(h.β)")
      for j = 1:N
        print(io, " $(-h.a[j])")
      end
      println(io)
    end
  else
    for v in vrep(rep)
      print(io, " $(Int(isray(v)))")
      for j = 1:N
        print(io, " $(v[j])")
      end
      println(io)
    end
  end
  print(io, "end")
end
