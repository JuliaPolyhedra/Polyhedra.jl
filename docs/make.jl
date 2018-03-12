using Documenter, Polyhedra

makedocs(
    format = :html,
    sitename = "Polyhedra",
    pages = [
        "Index" => "index.md",
        "Installation" => "installation.md",
        "Representation" => "representation.md",
        "Polyhedron" => "polyhedron.md",
        "Containment/Redundancy" => "redundancy.md",
        "Projection/Elimination" => "projection.md",
        "Utilities" => "utilities.md"
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [Polyhedra]
)

deploydocs(
    repo   = "github.com/JuliaPolyhedra/Polyhedra.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
