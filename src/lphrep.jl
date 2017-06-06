export LPHRepresentation

# No copy since I do not modify anything and a copy is done when building a polyhedron
mutable struct LPHRepresentation{N, T} <: HRepresentation{N, T}
    # lb <= Ax <= ub
    # l <= x <= u
    A::AbstractMatrix{T}
    l::AbstractVector{T}
    u::AbstractVector{T}
    colleqs::IntSet
    colgeqs::IntSet
    coleqs::IntSet
    lb::AbstractVector{T}
    ub::AbstractVector{T}
    rowleqs::IntSet
    rowgeqs::IntSet
    roweqs::IntSet

    function LPHRepresentation{N, T}(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) where {N, T}
        if length(l) != length(u) || size(A, 2) != length(l)
            error("The length of l and u must be equal to the number of rows of A")
        end
        if length(lb) != length(ub) || size(A, 1) != length(lb)
            error("The length of lb and ub must be equal to the number of columns of A")
        end
        if size(A, 2) != N
            error("Type dimension does not match the number of rows of A")
        end
        colleqs = IntSet()
        colgeqs = IntSet()
        coleqs = IntSet()
        for i in 1:N
            leq = u[i] < typemax(T)
            geq = l[i] > typemin(T)
            if leq && geq && myeq(l[i], u[i])
                push!(coleqs, i)
            else
                if leq
                    push!(colleqs, i)
                end
                if geq
                    push!(colgeqs, i)
                end
            end
        end
        rowleqs = IntSet()
        rowgeqs = IntSet()
        roweqs = IntSet()
        for i in 1:size(A, 1)
            leq = ub[i] < typemax(T)
            geq = lb[i] > typemin(T)
            if leq && geq && myeq(lb[i], ub[i])
                push!(roweqs, i)
            else
                if leq
                    push!(rowleqs, i)
                end
                if geq
                    push!(rowgeqs, i)
                end
            end
        end
        new{N, T}(A, l, u, colleqs, colgeqs, coleqs, lb, ub, rowleqs, rowgeqs, roweqs)
    end
end

decomposedfast(lp::LPHRepresentation) = false

LPHRepresentation{T <: Real}(A::AbstractMatrix{T}, l::AbstractVector{T}, u::AbstractVector{T}, lb::AbstractVector{T}, ub::AbstractVector{T}) = LPHRepresentation{size(A,2),T}(A, l, u, lb, ub)
function LPHRepresentation(A::AbstractMatrix, l::AbstractVector, u::AbstractVector, lb::AbstractVector, ub::AbstractVector)
    T = promote_type(eltype(A), eltype(l), eltype(u), eltype(lb), eltype(ub))
    LPHRepresentation{size(A,2),T}(AbstractMatrix{T}(A), AbstractVector{T}(l), AbstractVector{T}(u), AbstractVector{T}(lb), AbstractVector{T}(ub))
end

LPHRepresentation(rep::LPHRepresentation) = rep
LPHRepresentation{N,T}(rep::HRep{N,T}) = LPHRepresentation{N,T}(rep)

function (::Type{LPHRepresentation{N, T}}){N,T}(it::HRepIterator{N, T})
    A = Matrix{T}(length(it), N)
    lb = Vector{T}(length(it))
    ub = Vector{T}(length(it))
    MathProgBase.HighLevelInterface.warn_no_inf(T)
    l = fill(typemin(T), N)
    u = fill(typemax(T), N)
    for (i, h) in enumerate(it)
        A[i,:] = h.a
        ub[i] = h.β
        if islin(h)
            lb[i] = ub[i]
        else
            lb[i] = typemin(T)
        end
    end
    LPHRepresentation{N, T}(A, l, u, lb, ub)
end
function (::Type{LPHRepresentation{N, T}}){N,T}(eqs, ineqs)
    neq = length(eqs)
    nhrep = neq + length(ineqs)
    A = Matrix{T}(nhrep, N)
    lb = Vector{T}(nhrep)
    ub = Vector{T}(nhrep)
    MathProgBase.HighLevelInterface.warn_no_inf(T)
    l = fill(typemin(T), N)
    u = fill(typemax(T), N)
    for (i, h) in enumerate(eqs)
        A[i,:] = h.a
        lb[i] = h.β
        ub[i] = h.β
    end
    for (i, h) in enumerate(ineqs)
        A[neq+i,:] = h.a
        lb[neq+i] = typemin(T)
        ub[neq+i] = h.β
    end
    LPHRepresentation{N, T}(A, l, u, lb, ub)
end

Base.copy{N,T}(lp::LPHRepresentation{N,T}) = LPHRepresentation{N,T}(copy(A), copy(l), copy(u), copy(colleqs), copy(colgeqs), copy(coleqs), copy(lb), copy(ub), copy(rowleqs), copy(rowgeqs), copy(roweqs))

function checknext(lp::LPHRepresentation, state, allowed)
    colrow, i, lgeq = state
    lgeq += 1
    ok = false
    while colrow <= 2 && !ok
        if colrow == 1
            lgeqs = (lp.colleqs, lp.colgeqs, lp.coleqs)
        else
            lgeqs = (lp.rowleqs, lp.rowgeqs, lp.roweqs)
        end
        while i <= (colrow == 1 ? size(lp.A, 2) : size(lp.A, 1)) && !ok
            while lgeq <= 3 && !ok
                if allowed(lgeq) && i in lgeqs[lgeq]
                    ok = true
                else
                    lgeq += 1
                end
            end
            if !ok
                i += 1
                lgeq = 1
            end
        end
        if !ok
            colrow += 1
            i = 1
            lgeq = 1
        end
    end
    (colrow, i, lgeq)
end

neqs(lp::LPHRepresentation) = length(lp.coleqs) + length(lp.roweqs)
nineqs(lp::LPHRepresentation) = length(lp.colleqs) + length(lp.colgeqs) + length(lp.rowleqs) + length(lp.rowgeqs)

function gethrepaux{N, T}(lp::LPHRepresentation{N, T}, state)
    colrow, i, lgeq = state
    if colrow == 1
        a = spzeros(T, N)
        a[i] = lgeq == 2 ? -one(T) : one(T)
        β = lgeq == 2 ? -lp.l[i] : lp.u[i]
    elseif colrow == 2
        a = lgeq == 2 ? -lp.A[i,:] : lp.A[i,:]
        β = lgeq == 2 ? -lp.lb[i] : lp.ub[i]
    else
        error("The iterator is done")
    end
    lgeq == 3 ? HyperPlane(a, β) : HalfSpace(a, β)
end

starthrep(lp::LPHRepresentation) = checknext(lp, (1, 0, 3), (i) -> true)
donehrep(lp::LPHRepresentation, state) = state[1] > 2
function nexthrep(lp::LPHRepresentation, state)
    (gethrepaux(lp, state), checknext(lp, state, (i) -> true))
end

starteq(lp::LPHRepresentation) = checknext(lp, (1, 0, 3), (i) -> i == 3)
doneeq(lp::LPHRepresentation, state) = state[1] > 2
function nexteq(lp::LPHRepresentation, state)
    (gethrepaux(lp, state), checknext(lp, state, (i) -> i == 3))
end

startineq(lp::LPHRepresentation) = checknext(lp, (1, 0, 3), (i) -> i <= 2)
doneineq(lp::LPHRepresentation, state) = state[1] > 2
function nextineq(lp::LPHRepresentation, state)
    (gethrepaux(lp, state), checknext(lp, state, (i) -> i <= 2))
end
