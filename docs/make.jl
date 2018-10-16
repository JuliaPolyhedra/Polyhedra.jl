# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

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
        "Optimization" => "optimization.md",
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
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
