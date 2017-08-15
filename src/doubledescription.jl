export doubledescription
# Naive implementation of the double description Algorithm
# See JuliaPolyhedra/ConvexHull.jl for an efficient implementation

polytypefor{T <: Integer}(::Type{T}) = Rational{BigInt}
polytypefor(::Type{Float32}) = Float64
polytypefor{T}(::Type{T}) = T

function dualfullspace(h::HRepresentation{N, Tin}) where {N, Tin}
    Tout = polytypefor(Tin)
    SimpleVRepresentation(Matrix{Tout}(0, N), eye(Tout, N), IntSet(), IntSet(1:N))
end

function doubledescription(h::HRepresentation)
    # The redundancy of V-elements are measured using
    # the number of hyperplane they are in. If there are
    # redundant hyperplanes, it does not matter since no point
    # will be inside them but if there are duplicates it is a problem

    h = removeduplicates(h)
    v = dualfullspace(h)
    for hel in hreps(h)
        #v = removeduplicates(v ∩ hel)
        # removeduplicates is cheaper than removevredundancy since the latter
        # needs to go through all the hrep element
        v = removevredundancy(removeduplicates(v ∩ hel), h)
        #v = removevredundancy(v ∩ hel, h)
    end
    v
end

function doubledescription(v::VRepresentation{N, T}) where {N, T}
    lv = LiftedVRepresentation(v)
    # See #28
    if haspoints(v)
        R = -lv.R
    else
        R = [-one(T) zeros(T, 1, N); -lv.R]
    end
    vl = doubledescription(SimpleHRepresentation(R, zeros(T, size(R, 1)), lv.linset))
    # it has been cut by homogeneous hyperplane (i.e. containing the origin) only
    @assert npoints(vl) == 0
    LiftedHRepresentation(vl.R, vl.Rlinset)
end
