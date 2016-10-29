include("perm.jl")

using GLVisualize

window = glscreen()
view(visualize(GLNormalMesh(poly2)), window)
renderloop(window)
