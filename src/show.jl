import Base.show

_length(hrep::HRepresentation) = nhyperplanes(hrep) + nhalfspaces(hrep)
_length(vrep::VRepresentation) = nsympoints(vrep) + npoints(vrep) + nlines(vrep) + nrays(vrep)

function _print(io::IO, it::ElemIt{<:HRepElement})
    for h in it
        print(io, " $(h.β)")
        for j in eachindex(h.a)
            print(io, " $(-h.a[j])")
        end
        println(io)
    end
end

function _print(io::IO, it::ElemIt{<:VRepElement})
    for v in it
        print(io, " $(Int(ispoint(v)))")
        c = coord(v)
        for j = eachindex(c)
            print(io, " $(c[j])")
        end
        println(io)
    end
end

function Base.show(io::IO, rep::Representation{N,T}) where {N,T}
    if typeof(rep) <: HRepresentation
        print(io, "H")
        ls = IntSet(1:nhyperplanes(rep))
    else
        print(io, "V")
        ls = IntSet(1:nsympoints(rep)) ∪ IntSet(nsympoints(rep) + npoints(rep) + (1:nlines(rep)))
    end
    println(io, "-representation")

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
    println(io, " $(_length(rep)) $(N+1) $typename")
    if typeof(rep) <: HRepresentation
        _print.(io, hreps(rep))
    else
        _print.(io, vreps(rep))
    end
    print(io, "end")
end
