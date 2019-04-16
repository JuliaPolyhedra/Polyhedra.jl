using Documenter, Polyhedra

makedocs(
    sitename = "Polyhedra",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Installation" => "installation.md",
        "Representation" => "representation.md",
        "Polyhedron" => "polyhedron.md",
        "Plot" => "plot.md",
        "Containment/Redundancy" => "redundancy.md",
        "Projection/Elimination" => "projection.md",
        "Optimization" => "optimization.md",
        "Utilities" => "utilities.md"
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [Polyhedra]
)

deploydocs(
    repo   = "github.com/JuliaPolyhedra/Polyhedra.jl.git",
)
