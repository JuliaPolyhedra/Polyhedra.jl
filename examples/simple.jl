using Polyhedra
using CDDLib
using GLVisualize
using GeometryTypes

ext1 = SimpleVRepresentation([0 1 2;0 2 1; 1 0 2; 1 2 0; 2 0 1; 2 1 0])
ext2 = SimpleVRepresentation([0 0 0;0 0 1; 0 1 0; 1 0 0])
ext3 = SimpleVRepresentation([0 0 0;0 0 1; 0 1 0; 1 0 0], IntSet([1])) # BUG in CDD : 0x + 0y + 0z <= 1 -> useless
poly1 = CDDPolyhedron{3,Rational{BigInt}}(ext1)
poly2 = CDDPolyhedron{3,Rational{BigInt}}(ext2)
poly3 = CDDPolyhedron{3,Rational{BigInt}}(ext3)

window = glscreen()

#view(visualize(GeometryTypes.GLNormalMesh(poly1, RGBA(1f0,0f0,0f0,1f0))), window)
view(visualize(GLNormalMesh(poly1)), window)
#view(visualize(poly2), window)
#view(visualize(GeometryTypes.GLNormalMesh(poly3, RGBA(0f0,1f0,0f0,1f0))), window)
view(visualize(GLNormalMesh(poly3)), window)

include("axes.jl")

renderloop(window)
