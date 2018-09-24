# Show for elements
Base.show(io::IO, h::HyperPlane) = print(io, "HyperPlane($(h.a), $(h.β))")
Base.show(io::IO, h::HalfSpace) = print(io, "HalfSpace($(h.a), $(h.β))")
Base.show(io::IO, l::Line) = print(io, "Line($(l.a))")
Base.show(io::IO, r::Ray) = print(io, "Ray($(r.a))")

Base.summary(it::Polyhedra.AbstractRepIterator{T, ElemT}) where {T, ElemT} = "$(length(it))-element iterator of $ElemT"

# Inspired from Base.show_vector
function show_repit(io::IO, v::Polyhedra.AbstractRepIterator, print_prefix::Bool, start_str::String, end_str::String, join_str::String=",")
    if print_prefix
        print(io, Base.typeinfo_prefix(io, v))
    end
    io = IOContext(io, :typeinfo => eltype(v), :compact => get(io, :compact, true))
    limited = get(io, :limit, false)
    limited = true
    if limited && length(v) > 20
        Base.show_delim_array(io, v, start_str, join_str, " \u2026 " * end_str, false, 1, 20)
    else
        Base.show_delim_array(io, v, start_str, join_str, end_str, false)
    end
end

# Inspired from Base.showarray with repr=false and Base.print_matrix
function show_repit(io::IO, it::Polyhedra.AbstractRepIterator, ::MIME"text/plain")
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


Base.show(io::IO, v::Polyhedra.AbstractRepIterator) = show_repit(io, v, true, "[", "]")
Base.show(io::IO, mime::MIME"text/plain", it::Polyhedra.AbstractRepIterator) = show_repit(io, it, mime)

Base.summary(v::VRepresentation) = "V-representation $(typeof(v))"
Base.summary(h::HRepresentation) = "H-representation $(typeof(h))"
Base.summary(p::Polyhedron) = "Polyhedron $(typeof(p))"

show_reps(io::IO, args::Tuple, start_str::String, join_str::String, first::Bool) = first
function show_reps(io::IO, args::Tuple, start_str::String, join_str::String, first::Bool, rep, reps...)
    if !isempty(rep)
        if first
            print(io, start_str)
            first = false
        else
            print(io, join_str)
        end
        show_repit(io, rep, args...)
    end
    show_reps(io::IO, args, start_str, join_str, first, reps...)
end

show_vreps(io::IO, rep::HRepresentation, start_str::String, join_str::String, args...) = true
show_vreps(io::IO, rep::VRepresentation, start_str::String, join_str::String, args...) = show_reps(io, args, start_str, join_str, true, vreps(rep)...)
function show_vreps(io::IO, rep::Polyhedron, start_str::String, join_str::String, args...)
    if vrepiscomputed(rep)
        show_reps(io, args, start_str, join_str, true, vreps(rep)...)
    else
        true
    end
end
show_hreps(io::IO, rep::VRepresentation, start_str::String, join_str::String, args...) = true
show_hreps(io::IO, rep::HRepresentation, start_str::String, join_str::String, args...) = show_reps(io, args, start_str, join_str, true, hreps(rep)...)
function show_hreps(io::IO, rep::Polyhedron, start_str::String, join_str::String, args...)
    if hrepiscomputed(rep)
        show_reps(io, args, start_str, join_str, true, hreps(rep)...)
    else
        true
    end
end

_has_vrep(p::VRepresentation) = true
_has_vrep(p::HRepresentation) = false
_has_vrep(p::Polyhedron) = vrepiscomputed(p)

function Base.show(io::IO, rep::Rep)
    first = show_hreps(io, rep, "", " ∩ ", false, "", "", " ∩")
    if _has_vrep(rep) && !first
        print(io, " : ")
    end
    show_vreps(io, rep, "", " + ", false, "convexhull(", ")")
end
function Base.show(io::IO, mime::MIME"text/plain", rep::Rep)
    print(io, summary(rep))
    show_hreps(io, rep, ":\n", ",\n", mime)
    show_vreps(io, rep, ":\n", ",\n", mime)
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
#function Base.show(io::IO, rep::Representation{T}) where {T}
#    if typeof(rep) <: HRepresentation
#        print(io, "H")
#        ls = BitSet(1:nhyperplanes(rep))
#    else
#        print(io, "V")
#        ls = BitSet(1:nsympoints(rep)) ∪ BitSet(nsympoints(rep) + npoints(rep) + (1:nlines(rep)))
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
