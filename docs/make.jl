using Polyhedra
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = [
    "Convex hull and intersection",
    "Extended Formulation",
    "Minimal Robust Positively Invariant Set",
    "Convex hull of a set of points",
    "Projection of H-representation",
]

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example * ".jl")
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
end

# See https://juliadocs.github.io/Documenter.jl/v0.25/man/doctests/#Setup-Code
DocMeta.setdocmeta!(Polyhedra, :DocTestSetup, :(using Polyhedra); recursive=true)

makedocs(
    sitename = "Polyhedra",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Index" => "index.md",
        "Installation" => "installation.md",
        "Representation" => "representation.md",
        "Polyhedron" => "polyhedron.md",
        "Plot" => "plot.md",
        "Containment/Redundancy" => "redundancy.md",
        "Projection/Elimination" => "projection.md",
        "Optimization" => "optimization.md",
        "Utilities" => "utilities.md",
        "Internal" => "internal.md",
        "Examples" => map(EXAMPLES) do example
            example => "generated/$example.md"
        end
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [Polyhedra]
)

deploydocs(
    repo   = "github.com/JuliaPolyhedra/Polyhedra.jl.git",
    push_preview = true,
)
