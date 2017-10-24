using Documenter, Polyhedra

makedocs(
    format = :html,
    sitename = "Polyhedra",
    pages = [
        "Index" => "index.md",
        "Installation" => "installation.md",
        "Representation" => "representation.md",
        "Polyhedron" => "polyhedron.md"
    ]
)

deploydocs(
    repo   = "github.com/JuliaPolyhedra/Polyhedra.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
