export doubledescription
# Naive implementation of the double description Algorithm
# See JuliaPolyhedra/ConvexHull.jl for an efficient implementation

polytypefor{T <: Integer}(::Type{T}) = Rational{BigInt}
polytypefor(::Type{Float32}) = Float64
polytypefor{T}(::Type{T}) = T

function dualfullspace(h::HRepresentation{N, Tin}) where {N, Tin}
    Tout = polytypefor(Tin)
    MixedMatVRep(zeros(Tout, 1, N), eye(Tout, N), IntSet(), IntSet(1:N))
end

function doubledescription(h::HRepresentation)
    # The redundancy of V-elements are measured using
    # the number of hyperplane they are in. If there are
    # redundant hyperplanes, it does not matter since no point
    # will be inside them but if there are duplicates it is a problem

    h = removeduplicates(h)
    v = dualfullspace(h)
    checkvconsistency(v)
    for hp in hyperplanes(h)
        #v = removeduplicates(v ∩ hel)
        # removeduplicates is cheaper than removevredundancy since the latter
        # needs to go through all the hrep element
        v = removevredundancy(removeduplicates(v ∩ hp), h)
        #v = removevredundancy(v ∩ hel, h)
    end
    for hs in halfspaces(h)
        #v = removeduplicates(v ∩ hel)
        # removeduplicates is cheaper than removevredundancy since the latter
        # needs to go through all the hrep element
        v = removevredundancy(removeduplicates(v ∩ hs), h)
        #v = removevredundancy(v ∩ hel, h)
    end
    v
end

function doubledescription(v::VRepresentation{N, T}) where {N, T}
    checkvconsistency(v)
    lv = LiftedVRepresentation(v)
    R = -lv.R
    vl = doubledescription(MixedMatHRep(R, zeros(T, size(R, 1)), lv.linset))
    LiftedHRepresentation(vl.R, vl.Rlinset)
end
