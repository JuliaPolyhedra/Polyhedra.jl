using Polyhedra
using Base.Test
using CDD
using GLVisualize
using GeometryTypes
using Colors

ext1 = GeneratorDescription([0 1 2;0 2 1; 1 0 2; 1 2 0; 2 0 1; 2 1 0])
ext2 = GeneratorDescription([0 0 0;0 0 1; 0 1 0; 1 0 0])
ext3 = GeneratorDescription([0 0 0;0 0 1; 0 1 0; 1 0 0], IntSet([1])) # BUG in CDD : 0x + 0y + 0z <= 1 -> useless
poly1 = CDDPolyhedra(ext1)
poly2 = CDDPolyhedra(ext2)
poly3 = CDDPolyhedra(ext3)

window = glscreen()
view(visualize(GeometryTypes.GLNormalMesh(poly1, RGBA(1f0,0f0,0f0,1f0))), window)
#view(visualize(poly2), window)
view(visualize(GeometryTypes.GLNormalMesh(poly3, RGBA(0f0,1f0,0f0,1f0))), window)


# http://www.glvisualize.com/examples/meshes/
baselen = 0.1f0
dirlen = 2f0
# create an array of differently colored boxes in the direction of the 3 axes
rectangles = [
(HyperRectangle{3,Float32}(Vec3f0(-baselen/2), Vec3f0(baselen)), RGBA(0f0,0f0,0f0,1f0)),
(HyperRectangle{3,Float32}(Vec3f0(baselen/2, -baselen/2, -baselen/2), Vec3f0(dirlen, baselen, baselen)), RGBA(1f0,0f0,0f0,1f0)),
(HyperRectangle{3,Float32}(Vec3f0(-baselen/2, baselen/2, -baselen/2), Vec3f0(baselen, dirlen, baselen)), RGBA(0f0,1f0,0f0,1f0)),
(HyperRectangle{3,Float32}(Vec3f0(-baselen/2, -baselen/2, baselen/2), Vec3f0(baselen, baselen, dirlen)), RGBA(0f0,0f0,1f0,1f0))
]
# convert to an array of normal meshes
# note, that the constructor is a bit weird. GLNormalMesh takes a tuple of
# a geometry and a color. This means, the geometry will be converted to a GLNormalMesh
# and the color will be added afterwards, so the resulting type is a GLNormalColorMesh
meshes = map(GLNormalMesh, rectangles)
# merge them into one big mesh
# the resulting type is a GLNormalAttributeMesh, since we merged meshes with different
# attributes (colors). An array of the colors will be created and each vertex in the
# mesh will be asigned to one of the colors found there.
colored_mesh = merge(meshes)
view(visualize(colored_mesh), window)

renderloop(window)
