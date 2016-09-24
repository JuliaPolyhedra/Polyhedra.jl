using MathProgBase

import MathProgBase.LinearQuadraticModel, MathProgBase.loadproblem!, MathProgBase.optimize!, MathProgBase.status, MathProgBase.getobjval, MathProgBase.getsolution, MathProgBase.getunboundedray, MathProgBase.linprog

export AbstractPolyhedraModel

abstract AbstractPolyhedraModel <: MathProgBase.AbstractLinearQuadraticModel

loadproblem!(m::AbstractPolyhedraModel, hrep::HRep, c, sense) = error("loadproblem! not implemented")

include("lphrep.jl")

function loadproblem!(m::AbstractPolyhedraModel, A, l, u, c, lb, ub, sense)
  @show A, l, u, c, lb, ub
  loadproblem!(m, LPHRepresentation(A, l, u, lb, ub), c, sense)
end
function loadproblem!(m::MathProgBase.AbstractLinearQuadraticModel, hrep::HRep, c, sense)
  lp = LPHRepresentation(hrep)
  loadproblem!(p, lp.A, lp.l, lp.u, c, lp.lb, lp.ub, sense)
end

export SimpleVRepPolyhedraModel

type SimpleVRepPolyhedraModel <: AbstractPolyhedraModel
  vrep::Nullable{VRep}
  obj::Nullable{Vector}
  sense::Nullable{Symbol}

  objval::Nullable
  solution::Nullable{Vector}
  status::Symbol

  function SimpleVRepPolyhedraModel()
    new(nothing, nothing, nothing, nothing, nothing, :Undecided)
  end
end

function loadproblem!(lpm::SimpleVRepPolyhedraModel, vrep::VRep, obj, sense)
  if !(sense in [:Max, :Min])
    error("sense should be :Max or :Min")
  end
  if !myeqzero(obj)
    lpm.obj = copy(obj)
    lpm.sense = sense
  end
end

function myobjval(obj, v::VRepElement)
  if islin(r)
    if lpm.sense == :Min
      objval = min(mydot(lpm.obj, coord(r)), -mydot(lpm.obj, coord(r)))
    else
      objval = max(mydot(lpm.obj, coord(r)), -mydot(lpm.obj, coord(r)))
    end
  else
    objval = mydot(lpm.obj, coord(r))
  end
end

function optimize!(lpm::SimpleVRepPolyhedraModel)
  if isnull(lpm.vrep)
    error("Not problem loaded")
  end
  prob = get(lpm.vrep)
  N = fulldim(prob)
  T = eltype(prob)
  if !hasvreps(prob)
    lpm.status = :Infeasible
  else
    if lpm.sense in [:Max, :Min]
      better(a, b) = (lpm.sense == :Max ? a > b : a < b)
      mybetter(a, b) = (lpm.sense == :Max ? mygt(a, b) : mylt(a, b))
      lpm.status = :Infeasible
      for r in rays(prob)
        objval = myobjval(lpm.obj, r)
        if lpm.status != :Unbounded && mybetter(objval, zero(T))
          lpm.status = :Unbounded
          lpm.objval = lpm.sense == :Max ? typemax(T) : typemin(T)
          lpm.solution = r
        end
      end
      if status != :Unbounded
        for p in points(prob)
          objval = myobjval(lpm.obj, p)
          if lpm.status == :Undecided || (lpm.status == :Optimal && better(objval, get(lpm.objval)))
            lpm.status = :Optimal
            lpm.objval = objval
            lpm.solution = p
          end
        end
      end
    else
      lpm.status = :Optimal
      lpm.objval = zero(T)
      lpm.solution = spzeros(T,N)
    end
  end
end

function MathProgBase.status(lpm::SimpleVRepPolyhedraModel)
  lpm.status
end
function MathProgBase.getobjval(lpm::SimpleVRepPolyhedraModel)
  get(lpm.objval)
end
function MathProgBase.getsolution(lpm::SimpleVRepPolyhedraModel)
  copy(get(lpm.solution))
end
function MathProgBase.getunboundedray(lpm::SimpleVRepPolyhedraModel)
  copy(get(lpm.solution))
end

type LinprogSolution
    status
    objval
    sol
    attrs
end

function MathProgBase.linprog{N}(c::Vector, p::Rep{N}, solver::MathProgBase.AbstractMathProgSolver = defaultLPsolverfor(p))
  if N != length(c)
    println("length of objective does not match dimension of polyhedron")
  end
  loadproblem!(m, c, :Min)
  optimize!(m)
  stat = status(m)
  if stat == :Optimal
    return LinprogSolution(stat, getobjval(m), getsolution(m), Dict())
  elseif stat == :Unbounded
    attrs = Dict()
    try
      attrs[:unboundedray] = getunboundedray(m)
    catch
      warn("Problem is unbounded, but unbounded ray is unavailable; check that the proper solver options are set.")
    end
    return LinprogSolution(stat, nothing, [], attrs)
  else
    return LinprogSolution(stat, nothing, [], Dict())
  end
end

function Base.isempty{N,T}(p::Polyhedron{N,T})
  linprog(zeros(T, N), p).status == :Infeasible
end

function ishredundantaux(p::Polyhedron, a, b, strongly, cert, solver)
  sol = linprog(-a, p, solver)
  if sol.status == :Unbounded
    cert ?  (false, sol.attrs[:unboundedray], :UnboundedRay) : false
  elseif sol.status == :Optimal
    if mygt(sol.objval, b)
      cert ? (false, sol.sol, :ExteriorPoint) : false
    elseif mygeq(sol.objval, b)
      if strongly
        cert ? (false, sol.sol, :BoundaryPoint) : false
      else
        cert ? (true, sol.sol, :BoundaryPoint) : true
      end
    else
      cert ? (true, nothing, :NotApplicable) : true
    end
  end
end
function ishredundant(p::Rep, h::HRepElement; strongly=false, cert=false, solver = defaultLPsolverfor(p))
  if islin(h)
    sol = ishredundantaux(p, h.a, h.β, strongly, cert, solver)
    if !sol[1]
      sol
    else
      ishredundantaux(p, -h.a, -h.β, strongly, cert, solver)
    end
  else
    ishredundantaux(p, h.a, h.β, strongly, cert, solver)
  end
end
