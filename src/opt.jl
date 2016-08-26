using MathProgBase

import MathProgBase.LinearQuadraticModel, MathProgBase.loadproblem!, MathProgBase.optimize!, MathProgBase.status, MathProgBase.getobjval, MathProgBase.getsolution, MathProgBase.getunboundedray, MathProgBase.linprog

abstract AbstractPolyhedraModel{N, T} <: MathProgBase.AbstractLinearQuadraticModel

loadproblem!(m::AbstractPolyhedraModel, hrep::HRep, c, sense) = error("loadproblem! not implemented")

include("lphrep.jl")

function loadproblem!(m::AbstractPolyhedraModel, A, l, u, c, lb, ub, sense)
  loadproblem!(m, LPHRepresentation(A, l, u, lb, ub), c, sense)
end
function loadproblem!(m::MathProgBase.AbstractLinearQuadraticModel, lp::LPHRepresentation, c, sense)
  loadproblem!(p, lp.A, lp.l, lp.u, c, lp.lb, lp.ub, sense)
end
function loadproblem!(m::MathProgBase.AbstractLinearQuadraticModel, hrep::HRep, c, sense)
  loadproblem!(p, LPHRepresentation(hrep), c, sense)
end

export LPPolyhedron, SimpleVRepPolyhedraModel

abstract LPPolyhedron{N, T}

type SimpleVRepPolyhedraModel{N, T} <: AbstractPolyhedraModel{N, T}
  vrep::Nullable{VRep{N, T}}
  obj::Nullable{Vector{T}}
  sense::Nullable{Symbol}

  objval::Nullable{T}
  solution::Nullable{Vector{T}}
  status::Symbol

  function SimpleVRepPolyhedraModel()
    new(nothing, nothing, nothing, nothing, nothing, :Undecided)
  end
end

function loadproblem!{N,T}(lpm::SimpleVRepPolyhedraModel{N,T}, vrep::VRep{N,T}, obj::Vector{T}, sense)
  if !(sense in [:Max, :Min])
    error("sense should be :Max or :Min")
  end
  if !myeqzero(obj)
    lpm.obj = copy(obj)
    lpm.sense = sense
  end
end
loadproblem!{N,T}(lpm::SimpleVRepPolyhedraModel{N,T}, vrep::VRep{N,T}, obj::Vector, sense) = loadproblem!(lpm, vrep, Vector{T}(obj), sense)
loadproblem!{N,T}(lpm::SimpleVRepPolyhedraModel{N,T}, vrep::VRep{N}, obj::Vector, sense) = loadproblem!(lpm, VRepresentation{N,T}(vrep), Vector{T}(obj), sense)

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

function optimize!{N, T}(lpm::SimpleVRepPolyhedraModel{N, T})
  if !hasvreps(lpm.vrep)
    lpm.status = :Infeasible
  else
    if lpm.sense in [:Max, :Min]
      better(a, b) = (lpm.sense == :Max ? a > b : a < b)
      mybetter(a, b) = (lpm.sense == :Max ? mygt(a, b) : mylt(a, b))
      lpm.status = :Infeasible
      for r in rays(lpm.vrep)
        objval = myobjval(lpm.obj, r)
        if lpm.status != :Unbounded && mybetter(objval, zero(T))
          lpm.status = :Unbounded
          lpm.objval = lpm.sense == :Max ? typemax(T) : typemin(T)
          lpm.solution = r
        end
      end
      if status != :Unbounded
        for p in points(lpm.vrep)
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
