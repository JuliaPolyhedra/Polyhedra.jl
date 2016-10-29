using Colors
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
