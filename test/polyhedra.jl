using Polyhedra

include("config.jl")

include("jump.jl")
include("misc.jl")
include("decompose.jl")

polyhedratests = Dict(
    "jump" => jumptest,
    "misc" => misctest,
    "decompose" => decomposetest
)

@polytestset polyhedra true
