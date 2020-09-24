using Test
import Literate

const EXAMPLES_DIR = @__DIR__
const OUTPUT_DIR   = joinpath(@__DIR__, "generated")

const EXAMPLES = [
    "Extended Formulation.jl",
    "Minimal Robust Positively Invariant Set.jl"
]

@testset "test_examples.jl" begin
    @testset "$(example)" for example in EXAMPLES
        Literate.script(joinpath(EXAMPLES_DIR, example), OUTPUT_DIR)
        include(joinpath(OUTPUT_DIR, example))
    end
end
