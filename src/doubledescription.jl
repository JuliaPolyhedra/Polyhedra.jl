export doubledescription
# Naive implementation of the double description Algorithm
# See JuliaPolyhedra/ConvexHull.jl for an efficient implementation

polytypefor(::Type{T}) where {T <: Integer} = Rational{BigInt}
polytypefor(::Type{Float32}) = Float64
polytypefor(::Type{T}) where {T} = T

polyvectortype(a) = a
# TODO sparse halfspaces does not mean sparse points
polyvectortype(::Type{<:SparseVector{T}}) where T = Vector{T}

dualtype(RepT::Type{<:Representation}) = dualtype(RepT, polyvectortype(vectortype(RepT)))
function dualfullspace(rep::Representation, d::FullDim, ::Type{T}) where T
    dualfullspace(rep, d, T, polyvectortype(similar_type(vectortype(rep), d, T)))
end
function dualfullspace(rep::Representation{T}) where T
    dualfullspace(rep, FullDim(rep), polytypefor(T))
end

"""
    doubledescription(h::HRepresentation)

Computes the V-representation of the polyhedron represented by `h` using the Double-Description algorithm [1, 2].

    doubledescription(V::VRepresentation)

Computes the H-representation of the polyhedron represented by `v` using the Double-Description algorithm [1, 2].

[1] Motzkin, T. S., Raiffa, H., Thompson, G. L. and Thrall, R. M.
The double description method
*Contribution to the Theory of Games*, *Princeton University Press*, **1953**

[2] Fukuda, K. and Prodon, A.
Double description method revisited
*Combinatorics and computer science*, *Springer*, **1996**, 91-111
"""
function doubledescription end

function print_v_summary(v::VRep)
    print("$(npoints(v)) points, $(nrays(v)) rays and $(nlines(v)) lines")
end

function intersect_and_remove_redundancy(v, hs, h; verbose=0)
    if eltype(hs) <: HalfSpace
        str = "halfspace"
    else
        str = "hyperplane"
    end
    for (i, hel) in enumerate(hs)
        if verbose >= 1
            println("Intersecting $str $i/$(length(hs))")
        end
        v_int = v ∩ hel
        if verbose >= 3
            print("Removing duplicates: ")
            print_v_summary(v_int)
            println(".")
        end
        # removeduplicates is cheaper than removevredundancy since the latter
        # needs to go through all the hrep element
        v_dup = removeduplicates(v_int)
        if verbose >= 3
            print("Removing redundancy: ")
            print_v_summary(v_dup)
            println(".")
        end
        v = removevredundancy(v_dup, h)
        if verbose >= 2
            print("After intersection:  ")
            print_v_summary(v)
            println(".")
        end
    end
    return v
end

function doubledescription(h::HRepresentation; kws...)
    # The redundancy of V-elements are measured using
    # the number of hyperplane they are in. If there are
    # redundant hyperplanes, it does not matter since no point
    # will be inside them but if there are duplicates it is a problem

    # FIXME Note that removevredundancy, uses `h` which contains all hyperplanes and halfspaces
    # already taken into account but also all the other ones. We should check that this
    # is the right idea.
    h = removeduplicates(h)
    v = dualfullspace(h)
    checkvconsistency(v)
    v = intersect_and_remove_redundancy(v, hyperplanes(h), h; kws...)
    v = intersect_and_remove_redundancy(v, halfspaces(h), h; kws...)
    return v
end

function doubledescription(v::VRepresentation{T}; kws...) where {T}
    checkvconsistency(v)
    lv = convert(LiftedVRepresentation{T, Matrix{T}}, v)
    R = -lv.R
    vl = doubledescription(MixedMatHRep{T}(R, zeros(T, size(R, 1)), lv.linset); kws...)
    LiftedHRepresentation{T}(vl.R, vl.Rlinset)
end
