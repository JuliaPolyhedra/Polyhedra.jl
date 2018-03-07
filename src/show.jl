Base.summary(it::Polyhedra.AbstractRepIterator{N, T, ElemT}) where {N, T, ElemT} = "$(length(it))-element iterator of $ElemT"

# Inspired from Base.show_vector
function Base.show(io::IO, v::Polyhedra.AbstractRepIterator)
    compact, prefix = Base.array_eltype_show_how(v)
    limited = get(io, :limit, false)
    if compact && !haskey(io, :compact)
        io = IOContext(io, :compact => compact)
    end
    print(io, prefix)
    limited = true
    if limited && length(v) > 20
        Base.show_delim_array(io, v, "[", ",", " \u2026 ]", false, 1, 20)
    else
        Base.show_delim_array(io, v, "[", ",", "]", false)
    end
end

# Inspired from Base.showarray with repr=false and Base.print_matrix
function Base.show(io::IO, ::MIME"text/plain", it::Polyhedra.AbstractRepIterator)
    print(io, summary(it))
    if !isempty(it)
        if !get(io, :limit, false)
            screenheight = screenwidth = typemax(Int)
        else
            sz = displaysize(io)
            # We use -4 as this is what is used in Base.print_matrix.
            # This probably accounts for the summary and the vertical dots
            screenheight, screenwidth = sz[1] - 4, sz[2]
        end
        if screenheight < length(it)
            Base.show_delim_array(io, it, ":\n ", "\n", "\n  \u22ee", false, 1, screenheight)
        else
            Base.show_delim_array(io, it, ":\n ", "\n", "", false)
        end
    end
end

#_length(hrep::HRepresentation) = nhyperplanes(hrep) + nhalfspaces(hrep)
#_length(vrep::VRepresentation) = nsympoints(vrep) + npoints(vrep) + nlines(vrep) + nrays(vrep)
#
#function _print(io::IO, it::ElemIt{<:HRepElement})
#    for h in it
#        print(io, " $(h.β)")
#        for j in eachindex(h.a)
#            print(io, " $(-h.a[j])")
#        end
#        println(io)
#    end
#end
#
#function _print(io::IO, it::ElemIt{<:VRepElement})
#    for v in it
#        print(io, " $(Int(ispoint(v)))")
#        c = coord(v)
#        for j = eachindex(c)
#            print(io, " $(c[j])")
#        end
#        println(io)
#    end
#end
#
#function Base.show(io::IO, rep::Representation{N,T}) where {N,T}
#    if typeof(rep) <: HRepresentation
#        print(io, "H")
#        ls = IntSet(1:nhyperplanes(rep))
#    else
#        print(io, "V")
#        ls = IntSet(1:nsympoints(rep)) ∪ IntSet(nsympoints(rep) + npoints(rep) + (1:nlines(rep)))
#    end
#    println(io, "-representation")
#
#    if !isempty(ls)
#        print(io, "linearity $(length(ls))");
#        for i in ls
#            print(io, " $i")
#        end
#        println(io)
#    end
#
#    println(io, "begin")
#    if T <: AbstractFloat
#        typename = "real"
#    elseif T <: Integer
#        typename = "integer"
#    else
#        typename = "rational"
#    end
#    println(io, " $(_length(rep)) $(N+1) $typename")
#    if typeof(rep) <: HRepresentation
#        _print.(io, hreps(rep))
#    else
#        _print.(io, vreps(rep))
#    end
#    print(io, "end")
#end
