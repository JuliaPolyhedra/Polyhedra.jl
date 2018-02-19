@test default_library(Polyhedra.FullDim{2}(), Float64) isa SimplePolyhedraLibrary
@test default_library(Polyhedra.FullDim{1}(), Float64) isa IntervalLibrary
