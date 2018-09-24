using JuMP
import MathProgBase

struct MockSolver <: MathProgBase.AbstractMathProgSolver end
struct MockModel <: AbstractPolyhedraModel end
Polyhedra.PolyhedraModel(::MockSolver) = MockModel()
MathProgBase.LinearQuadraticModel(s::MockSolver) = PolyhedraToLPQPBridge(Polyhedra.PolyhedraModel(s))
function MathProgBase.loadproblem!(::MockModel, lp, obj, sense)
    @test obj == [1.0, 1.0, 1.0, -2.0]
    @test sense == :Max
end
function MathProgBase.optimize!(::MockModel) end
MathProgBase.status(::MockModel) = :Optimal
MathProgBase.getsolution(::MockModel) = [0.0, 1.0, 0.0, 1.0]
MathProgBase.getobjval(::MockModel) = -1.0
MathProgBase.getconstrsolution(::MockModel) = [1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0, -1.0]
MathProgBase.getconstrduals(::MockModel) = [-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, -0.0]

function polyhedramodeltest(s::MathProgBase.AbstractMathProgSolver)
    m = Model(solver=s)
    n = 3
    @variable(m, 0 <= x[1:n] <= 1)
    @variable(m, y == 1)
    @constraint(m, sum(x[i] for i=1:n) == 1)
    @constraint(m, x[1] - x[2] <= 1)
    @constraint(m, x[2] - x[3] >= 1)
    @objective(m, Max, sum(x) - 2y)

    lp = hrep(m)
    coloffset, rowoffset = Polyhedra.computeoffsets(lp)
    @test coloffset == [[3, -4], [5, -6], [7, -8], [1]]
    @test rowoffset == [[2], [9], [-10]]

    status = solve(m)
    @test status == :Optimal
    @test MathProgBase.getsolution(m.internalModel) ≈ [0.0, 1.0, 0.0, 1.0]
    @test MathProgBase.getobjval(m.internalModel) ≈ -1.0
    @test MathProgBase.getconstrsolution(m.internalModel) ≈ [1.0, -1.0, 1.0]
    @test MathProgBase.getconstrduals(m.internalModel) ≈ [1.0, 0.0, 0.0]
end

@testset "PolyhedraToLPQPBridge" begin
    polyhedramodeltest(MockSolver())
end
