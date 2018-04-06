using Polyhedra

include("jump.jl")
include("misc.jl")

polyhedratests = Dict("jump" => jumptest,
                      "misc" => misctest)

@polytestset polyhedra true
