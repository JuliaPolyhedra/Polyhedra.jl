var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "#Polyhedra-–-Manipulation-of-Polyhedra-in-Julia-1",
    "page": "Index",
    "title": "Polyhedra –- Manipulation of Polyhedra in Julia",
    "category": "section",
    "text": "Polyhedra is a package for polyhedra manipulations in Julia. It provides an unified interface for Polyhedra Manipulation Libraries such as CDDLib.jl, LRSLib.jl and QHull.Polyhedra can either be represented by a set of linear inequalities or by vertices and rays. In the first case, the points of the polyhedron are the points that satisfies all the inequalities and in the second case they are the points that can be expressed as a convex combination of the vertices plus a conic combination of the rays. The manipulations that Polyhedra can perform includeProjection: Projection of a polyhedron on a lower dimensional space, e.g. Fourier-Motzkin elimination.\nChanging the Representation\nVertex enumeration problem: Computing the extremal vertices and rays from an inequality representation\nConvex hull problem: Computing a set of linear inequalities describing the polyhedron from a vertex/ray representation\nRemoval of redundant inequalities or redundant vertices/rays.\nPlotting of 2D polyhedra using Plots\nDecomposition of 3D polyhedra into vertices and triangular faces, enabling easy visualization of 3D polyhedra using DrakeVisualizer or GLVisualize.Depending on the library, these manipulation can either be in floating point or exact rational arithmetic.Each operation has a default fallback implementation which is used in case the library does not support it. Polyhedra also includes a default library which does not implement anything, hence using every fallback.Polyhedra remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "#Installing-Polyhedra-1",
    "page": "Index",
    "title": "Installing Polyhedra",
    "category": "section",
    "text": "If you are familiar with Julia you can get started quickly by using the package manager to install Polyhedrajulia> Pkg.add(\"Polyhedra\")And a Polyhedra Manipulation Library, e.g.julia> Pkg.add(\"CDDLib\")"
},

{
    "location": "#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\", \"representation.md\", \"polyhedron.md\", \"plot.md\", \"redundancy.md\", \"projection.md\", \"optimization.md\", \"utilities.md\"]\nDepth = 2"
},

{
    "location": "installation/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "This section shows how to install Julia, Polyhedra and a Polyhedra Manipulation Library of your choice."
},

{
    "location": "installation/#Getting-Julia-1",
    "page": "Installation",
    "title": "Getting Julia",
    "category": "section",
    "text": "The first step is to install Julia. Polyhedra supports Julia v0.4 to Julia v1.0 but the latest version only supports Julia v0.7 and v1.0. Download links and more detailed instructions are available on the Julia website."
},

{
    "location": "installation/#Getting-Polyhedra-1",
    "page": "Installation",
    "title": "Getting Polyhedra",
    "category": "section",
    "text": "Open a Julia console (e.g. enter julia at the command line) and write] add PolyhedraTo start using Polyhedra, you can now just writejulia> using PolyhedraPolyhedra includes a default library supporting every operation but external library can also be used. See the next section on installing a library."
},

{
    "location": "installation/#Getting-Libraries-1",
    "page": "Installation",
    "title": "Getting Libraries",
    "category": "section",
    "text": "Many C libraries are are available for manipulating Polyhedra. Some of them works with floating point arithmetic and some of them can do the computation exactly using rational arithmetic and multiple precision libraries such as GMP. Julia also natively support Rational arithmetic using multiple precision libraries and of course floating point arithmetic. That makes the use of both arithmetic very easy and transparent. A list of Polyhedra Manipulation Libraries is available in the JuliaPolyhera website."
},

{
    "location": "representation/#",
    "page": "Representation",
    "title": "Representation",
    "category": "page",
    "text": ""
},

{
    "location": "representation/#Polyhedra.HRepresentation",
    "page": "Representation",
    "title": "Polyhedra.HRepresentation",
    "category": "type",
    "text": "HRepresentation{T<:Real}\n\nSupertype for H-representations with coefficient type T.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.VRepresentation",
    "page": "Representation",
    "title": "Polyhedra.VRepresentation",
    "category": "type",
    "text": "VRepresentation{T<:Real}\n\nSupertype for V-representations coefficient type T.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.Representation",
    "page": "Representation",
    "title": "Polyhedra.Representation",
    "category": "type",
    "text": "Representation{T<:Real}\n\nSupertype for H-(or V-)representations with coefficient type T.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.fulldim",
    "page": "Representation",
    "title": "Polyhedra.fulldim",
    "category": "function",
    "text": "fulldim(rep::Rep)::Int\n\nReturns the dimension of the space in which polyhedron, representation, element or vector is defined. That is, a straight line in a 3D space has fulldim 3 even if its dimension is 1.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.FullDim",
    "page": "Representation",
    "title": "Polyhedra.FullDim",
    "category": "constant",
    "text": "FullDim(p)::FullDim\n\nSimilar to fulldim but used for type stability with the vector type. If the vector type is StaticArrays.SVector then it returns a StaticArrays.Size.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.coefficient_type",
    "page": "Representation",
    "title": "Polyhedra.coefficient_type",
    "category": "function",
    "text": "coefficient_type(rep::Rep)\n\nReturns the type of the coefficients used in the representation of rep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Representation-1",
    "page": "Representation",
    "title": "Representation",
    "category": "section",
    "text": "Polyhedra can be described in 2 different ways.H-representation: As the intersection of finitely many halfspaces given by its facets.\nV-representation: As the convex hull of its vertices + the conic hull of its rays where \'+\' is the Minkowski sum.HRepresentation\nVRepresentation\nRepresentationIn Polyhedra.jl, those representations are given the respective abstract types HRepresentation and VRepresentation which are themself subtypes of Representation.These functions can be called on both H-representation and V-representationfulldim\nPolyhedra.FullDim\nPolyhedra.coefficient_type"
},

{
    "location": "representation/#Polyhedra.HalfSpace",
    "page": "Representation",
    "title": "Polyhedra.HalfSpace",
    "category": "type",
    "text": "struct HalfSpace{T, AT} <: HRepElement{T, AT}\n    a::AT\n    β::T\nend\n\nAn halfspace defined by the set of points x such that langle a x rangle le beta.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.HyperPlane",
    "page": "Representation",
    "title": "Polyhedra.HyperPlane",
    "category": "type",
    "text": "struct HyperPlane{T, AT} <: HRepElement{T, AT}\n    a::AT\n    β::T\nend\n\nAn hyperplane defined by the set of points x such that langle a x rangle = beta.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hrep",
    "page": "Representation",
    "title": "Polyhedra.hrep",
    "category": "function",
    "text": "hrep(p::Polyhedron)\n\nReturns an H-representation for the polyhedron p.\n\n\n\n\n\nhrep(hyperplanes::HyperPlaneIt; d::FullDim)\n\nCreates an affine space of full dimension d from the list of hyperplanes hyperplanes.\n\nExamples\n\nhrep([HyperPlane([0, 1, 0], 1), HyperPlane([0, 0, 1], 0)])\n\ncreates the 1-dimensional affine subspace containing all the points (x_1 0 0), i.e. the x_1-axis.\n\nhrep([HyperPlane([1, 1], 1), HyperPlane([1, 0], 0)])\n\ncreates the 0-dimensional affine subspace only containing the point (0 1).\n\n\n\n\n\nhrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt; d::FullDim)\n\nCreates an H-representation for the polyhedron of full dimension d equal to the intersection of the hyperplanes hyperplanes and halfspaces halfspaces.\n\nExamples\n\nFor instance, the simplex\n\nbeginalign*\n  x_1 + x_2 = 1 \n  x_1 geq 0 \n  x_2 geq 0\nendalign*\n\ncan be created as follows:\n\nhrep([HalfSpace([-1, 0], 0)], [HyperPlane([1, 1], 1), HalfSpace([0, -1], 0)])\n\n\n\n\n\nhrep(halfspaces::HalfSpaceIt; d::FullDim)\n\nCreates an H-representation for the polyhedron of full dimension d equal to the intersection of the halfspaces halfspaces.\n\nExamples\n\nFor instance, the polytope\n\nbeginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*\n\ncan be created as follows:\n\nhrep([HalfSpace([1, 1], 1), HalfSpace([1, -1], 0), HalfSpace([-1, 0], 0)])\n\n\n\n\n\nhrep(model::JuMP.Model)\n\nBuilds an H-representation from the feasibility set of the JuMP model model. Note that if non-linear constraint are present in the model, they are ignored.\n\n\n\n\n\nhrep(A::AbstractMatrix, b::AbstractVector, linset::BitSet=BitSet())\n\nCreates an H-representation for the polyhedron defined by the inequalities langle A_i x rangle = b_i if i in linset and langle A_i x rangle le b_i otherwise where A_i is the ith row of A, i.e. A[i,:] and b_i is b[i].\n\n\n\n\n\n"
},

{
    "location": "representation/#H-representation-1",
    "page": "Representation",
    "title": "H-representation",
    "category": "section",
    "text": "The fundamental element of an H-representation is the halfspaceHalfSpaceAn H-representation can be created as the intersection of several halfspaces. For instance, the polytopebeginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*can be created as follows:HalfSpace([1, 1], 1) ∩ HalfSpace([1, -1], 0) ∩ HalfSpace([-1, 0], 0)Even if HalfSpaces are enough to describe any polyhedron, it is sometimes important to represent the fact that the polyhedron is contained in an affine subspace. For instance, the simplexbeginalign*\n  x_1 + x_2 = 1 \n  x_1 geq 0 \n  x_2 geq 0\nendalign*is 1-dimensional even if it is defined in a 2-dimensional space.The fundamental element of an affine subspace is the hyperplaneHyperPlaneAn affine subspace can be created as the intersection of several hyperplanes. For instanceHyperPlane([1, 1], 1) ∩ HyperPlane([1, 0], 0)represents the 0-dimensional affine subspace only containing the point (0 1).To represent a polyhedron that is not full-dimensional, hyperplanes and halfspaces can be mixed in any order. For instance, the simplex defined above can be obtained as follows:HalfSpace([-1, 0], 0) ∩ HyperPlane([1, 1], 1) ∩ HalfSpace([0, -1], 0)In addition to being created incrementally with intersections, an H-representation can also be created using the hrep functionhrep"
},

{
    "location": "representation/#Polyhedra.halfspaces",
    "page": "Representation",
    "title": "Polyhedra.halfspaces",
    "category": "function",
    "text": "halfspaces(hrep::HRep)\n\nReturns an iterator over the halfspaces of the H-representation hrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.nhalfspaces",
    "category": "function",
    "text": "nhalfspaces(hrep::HRep)\n\nReturns the number of halfspaces of the H-representation hrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hashalfspaces",
    "page": "Representation",
    "title": "Polyhedra.hashalfspaces",
    "category": "function",
    "text": "hashalfspaces(hrep::HRep)\n\nReturns whether the H-representation hrep has any halfspace.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hyperplanes",
    "page": "Representation",
    "title": "Polyhedra.hyperplanes",
    "category": "function",
    "text": "hyperplanes(hrep::HRep)\n\nReturns an iterator over the hyperplanes of the H-representation hrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nhyperplanes",
    "page": "Representation",
    "title": "Polyhedra.nhyperplanes",
    "category": "function",
    "text": "nhyperplanes(hrep::HRep)\n\nReturns the number of hyperplanes of the H-representation hrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hashyperplanes",
    "page": "Representation",
    "title": "Polyhedra.hashyperplanes",
    "category": "function",
    "text": "hashyperplanes(hrep::HRep)\n\nReturns whether the H-representation hrep has any hyperplane.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.allhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.allhalfspaces",
    "category": "function",
    "text": "allhalfspaces(hrep::HRep)\n\nReturns an iterator over the halfspaces and hyperplanes in the H-representation hrep splitting hyperplanes in two halfspaces.\n\nExamples\n\nhrep = HyperPlane([1, 0], 1) ∩ HalfSpace([0, 1], 1)\ncollect(allhalfspaces(hrep)) # Returns [HalfSpace([1, 0]), HalfSpace([-1, 0]), HalfSpace([0, 1])]\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nallhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.nallhalfspaces",
    "category": "function",
    "text": "nallhalfspaces(hrep::HRep)\n\nReturns the number of halfspaces plus twice the number of hyperplanes in the H-representation hrep, i.e. length(allhalfspaces(hrep))\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hasallhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.hasallhalfspaces",
    "category": "function",
    "text": "hasallhalfspaces(hrep::HRep)\n\nReturns whether the H-representation hrep contains any halfspace or hyperplane.\n\n\n\n\n\n"
},

{
    "location": "representation/#Interface-1",
    "page": "Representation",
    "title": "Interface",
    "category": "section",
    "text": "An H-representation is represented as an intersection halfspaces and hyperplanes. The halfspaces can be obtained with halfspaces and the hyperplanes with hyperplanes. As an hyperplane langle a x rangle = beta is the intersection of the two halfspaces langle a x rangle le beta and langle a x rangle ge beta, even if the H-representation contains hyperplanes, a list of halfspaces whose intersection is the polyhedron can be obtained with allhalfspaces, which has nhalfspaces(p) + 2nhyperplanes(p) elements for an H-representation p since each hyperplane is split in two halfspaces.halfspaces\nnhalfspaces\nhashalfspaces\nhyperplanes\nnhyperplanes\nhashyperplanes\nallhalfspaces\nnallhalfspaces\nhasallhalfspaces"
},

{
    "location": "representation/#Polyhedra.Ray",
    "page": "Representation",
    "title": "Polyhedra.Ray",
    "category": "type",
    "text": "struct Ray{T, AT <: AbstractVector{T}}\n    a::AT\nend\n\nThe conic hull of a, i.e. the set of points λa where λ is any nonnegative real number.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.Line",
    "page": "Representation",
    "title": "Polyhedra.Line",
    "category": "type",
    "text": "struct Line{T, AT <: AbstractVector{T}}\n    a::AT\nend\n\nThe conic hull of a and -a, i.e. the set of points λa where λ is any real number.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.vrep",
    "page": "Representation",
    "title": "Polyhedra.vrep",
    "category": "function",
    "text": "vrep(p::Polyhedron)\n\nReturns a V-representation for the polyhedron p.\n\n\n\n\n\nvrep(lines::LineIt; d::FullDim)\n\nCreates an affine space of full dimension d from the list of lines lines.\n\nExamples\n\nvrep([Line([1, 0, 0]), Line([0, 1, 0])])\n\ncreates the 2-dimensional affine subspace containing all the points (x_1 x_2 0), i.e. the x_1x_2-plane.\n\n\n\n\n\nvrep(points::PointIt; d::FullDim)\n\nCreates a V-representation for the polytope of full dimension d equal to the convex hull of the points points.\n\nExamples\n\nThe convex hull of (0 0), (0 1) and (12 12) can be created as follows using exact arithmetic\n\nvrep([[0, 0], [0, 1], [1//2, 1//2]])\n\nor as follows using floating point arithmetic\n\nvrep([[0, 0], [0, 1], [1/2, 1/2]])\n\n\n\n\n\nvrep(lines::LineIt, rays::RayIt; d::FullDim)\n\nCreates a V-representation for the polyhedral cone of full dimension d equal to the conic hull of the lines lines and rays rays.\n\nExamples\n\nvrep([Line([0, 1])], [Ray([1, 0])])\n\ncreates a V-representation for the halfspace x_1 ge 0.\n\n\n\n\n\nvrep(rays::RayIt)\n\nCreates a V-representation for the polyhedral cone of full dimension d equal to the conic hull of the rays rays.\n\nExamples\n\nvrep([Ray([1, 0]), Ray([0, 1])])\n\ncreates a V-representation for positive orthant.\n\n\n\n\n\nvrep(points::PointIt, lines::LineIt, rays::RayIt)\n\nCreates a V-representation for the polyhedron of full dimension d equal to the minkowski sum of the convex hull of points with the conic hull of lines and rays.\n\n\n\n\n\nvrep(V::AbstractMatrix, R::AbstractMatrix, Rlinset::BitSet=BitSet())\n\nCreates a V-representation for the polyhedron defined by the points V_i, lines R_i if i in Rlinset and rays R_i otherwise where V_i (resp. R_i) is the ith row of V (resp. R), i.e. V[i,:] (resp. R[i,:]).\n\n\n\n\n\n"
},

{
    "location": "representation/#V-representation-1",
    "page": "Representation",
    "title": "V-representation",
    "category": "section",
    "text": "The fundamental elements of an V-representation are the points (represented by AbstractVectors and and raysRayA V-representation can be created as the minkowski sum between a convex hull of points and a conic hull of rays. For instance, the positive orthant without the simplex defined in the H-representation section can be created as follows:convexhull([1, 0], [0, 1]) + conichull([1, 0], [0, 1])The V-representation represents the polyhedron as a minkowski sum of a polytope and a polyhedral cone. The polytope is represented using a P-representation : a convex hull of points. The polyhedral cone is represented using an R-representation : a conic hull of rays.Even if rays are enough to describe any polyhedral cone, it is sometimes important to represent the fact that the polyhedron contains an affine subspace. For instance, the polyhedron created withconvexhull([1, 0], [0, 1]) + conichull([1, 1], [-1, -1])contains the line [1, 1].The fundamental element of an affine subspace is the lineLineAn affine subspace can be created as the conic hull/minkownski sum of several lines. For instanceconichull(Line([1, 0]), Line([0, 1]))represents the full space.In addition to being created incrementally with convex hull and minkowsky addition, a V-representation can also be created using the vrep functionvrep"
},

{
    "location": "representation/#Polyhedra.points",
    "page": "Representation",
    "title": "Polyhedra.points",
    "category": "function",
    "text": "points(vrep::VRep)\n\nReturns an iterator over the points of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.npoints",
    "page": "Representation",
    "title": "Polyhedra.npoints",
    "category": "function",
    "text": "npoints(vrep::VRep)\n\nReturns the number of points of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.haspoints",
    "page": "Representation",
    "title": "Polyhedra.haspoints",
    "category": "function",
    "text": "haspoints(vrep::VRep)\n\nReturns whether the V-representation vrep has any point.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.rays",
    "page": "Representation",
    "title": "Polyhedra.rays",
    "category": "function",
    "text": "rays(vrep::VRep)\n\nReturns an iterator over the rays of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nrays",
    "page": "Representation",
    "title": "Polyhedra.nrays",
    "category": "function",
    "text": "nrays(vrep::VRep)\n\nReturns the number of rays of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hasrays",
    "page": "Representation",
    "title": "Polyhedra.hasrays",
    "category": "function",
    "text": "hasrays(vrep::VRep)\n\nReturns whether the V-representation vrep has any ray.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.lines",
    "page": "Representation",
    "title": "Polyhedra.lines",
    "category": "function",
    "text": "lines(vrep::VRep)\n\nReturns an iterator over the lines of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nlines",
    "page": "Representation",
    "title": "Polyhedra.nlines",
    "category": "function",
    "text": "nlines(vrep::VRep)\n\nReturns the number of lines of the V-representation vrep.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.haslines",
    "page": "Representation",
    "title": "Polyhedra.haslines",
    "category": "function",
    "text": "haslines(vrep::VRep)\n\nReturns whether the V-representation vrep has any line.\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.allrays",
    "page": "Representation",
    "title": "Polyhedra.allrays",
    "category": "function",
    "text": "allrays(vrep::VRep)\n\nReturns an iterator over the rays and lines in the V-representation vrep splitting lines in two rays.\n\nExamples\n\nvrep = Line([1, 0]) + Ray([0, 1])\ncollect(allrays(vrep)) # Returns [Ray([1, 0]), Ray([-1, 0]), Ray([0, 1])]\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.nallrays",
    "page": "Representation",
    "title": "Polyhedra.nallrays",
    "category": "function",
    "text": "nallrays(vrep::VRep)\n\nReturns the number of rays plus twice the number of lines in the V-representation vrep, i.e. length(allrays(vrep))\n\n\n\n\n\n"
},

{
    "location": "representation/#Polyhedra.hasallrays",
    "page": "Representation",
    "title": "Polyhedra.hasallrays",
    "category": "function",
    "text": "hasallrays(vrep::VRep)\n\nReturns whether the V-representation vrep contains any ray or line.\n\n\n\n\n\n"
},

{
    "location": "representation/#Interface-2",
    "page": "Representation",
    "title": "Interface",
    "category": "section",
    "text": "A P-representation is represented as a convex hull points. The points can be obtained with points.points\nnpoints\nhaspointsAn R-representation is represented as a conic hull of lines and rays. The rays can be obtained with rays and the lines with lines. As a line r is the conic hull of of the two rays r and -r, even if the R-representation contains lines, a list of rays whose conic hull is the polyhedral cone can be obtained with allrays, which has nrays(R) + 2nlines(R) elements for an R-representation R since each line is split in two rays.rays\nnrays\nhasrays\nlines\nnlines\nhaslines\nallrays\nnallrays\nhasallrays"
},

{
    "location": "polyhedron/#",
    "page": "Polyhedron",
    "title": "Polyhedron",
    "category": "page",
    "text": ""
},

{
    "location": "polyhedron/#Polyhedra.doubledescription",
    "page": "Polyhedron",
    "title": "Polyhedra.doubledescription",
    "category": "function",
    "text": "doubledescription(h::HRepresentation)\n\nComputes the V-representation of the polyhedron represented by h using the Double-Description algorithm [1, 2].\n\ndoubledescription(V::VRepresentation)\n\nComputes the H-representation of the polyhedron represented by v using the Double-Description algorithm [1, 2].\n\n[1] Motzkin, T. S., Raiffa, H., Thompson, G. L. and Thrall, R. M. The double description method Contribution to the Theory of Games, Princeton University Press, 1953\n\n[2] Fukuda, K. and Prodon, A. Double description method revisited Combinatorics and computer science, Springer, 1996, 91-111\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.polyhedron",
    "page": "Polyhedron",
    "title": "Polyhedra.polyhedron",
    "category": "function",
    "text": "polyhedron(rep::Representation{T})\n\nCreates a polyhedron from the representation rep using the default library included in the Polyhedra package.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedron-1",
    "page": "Polyhedron",
    "title": "Polyhedron",
    "category": "section",
    "text": "As seen in the previous section, a polyhedron can be described in 2 ways: either using the H-representation (intersection of halfspaces) or the V-representation (convex hull of points and rays). The problem of computing the H-representation from the V-representation (or vice versa) is called the representation conversion problem. It can be solved by the Double-Description methoddoubledescriptionHowever, other methods exist such as the reverse search implemented by LRS and the quick hull algorithm implemented by qhull.This motivates the creation of a type representing polyhedra, transparently handling the conversion from H-representation to V-representation when needed for some operation. Just like the abstract type AbstractArray{N,T} represents an N-dimensional array with elements of type T, the abstract type Polyhedron{N,T} represents an N-dimensional polyhedron with elements of coefficient type T.There is typically one concrete subtype of Polyhedron by library. For instance, the CDD library defines CDDLib.Polyhedron and the LRS library defines LRSLib.Polyhedron. It must be said that the type T is not necessarily how the elements are stored internally by the library but the polyhedron will behave just like it is stored that way. For instance, when retreiving an H-(or V-)representation, the representation will be of type T. Therefore using Int for T may result in InexactError. For this reason, by default, the type T chosen is not a subtype of Integer.A polyhedron can be created from a representation and a library using the polyhedron function.polyhedronTo illustrate the usage of the polyhedron function, consider the following representations:hr = HalfSpace([1, 1], 1) ∩ HalfSpace([1, -1], 0) ∩ HalfSpace([-1, 0], 0)\nvre = convexhull([0, 0], [0, 1], [1//2, 1//2])\nvrf = convexhull([0, 0], [0, 1], [1/2, 1/2])One can use the CDD library, to create an instance of a concrete subtype of Polyhedron as follows:julia> using CDDLib\njulia> polyf = polyhedron(hr, CDDLib.Library())\njulia> typeof(polyhf)\nCDDLib.CDDLib.Polyhedron{2,Float64}We see that the library has choosen to deal with floating point arithmetic. This decision does not depend on the type of hr but only on the instance of CDDLib.Library given. CDDLib.Library creates CDDLib.Polyhedron of type either Float64 or Rational{BigInt}. One can choose the first one using CDDLib.Library(:float) and the second one using CDDLib.Library(:exact), by default it is :float.julia> poly = polyhedron(hr, CDDLib.Library(:exact))\njulia> typeof(poly)\nCDDLib.Polyhedron{2,Rational{BigInt}}The first polyhedron polyf can also be created from its V-representation using either of the 4 following lines:julia> polyf = polyhedron(vrf, CDDLib.Library(:float))\njulia> polyf = polyhedron(vrf, CDDLib.Library())\njulia> polyf = polyhedron(vre,  CDDLib.Library(:float))\njulia> polyf = polyhedron(vre,  CDDLib.Library())and poly using either of those lines:julia> poly = polyhedron(vrf, CDDLib.Library(:exact))\njulia> poly = polyhedron(vre, CDDLib.Library(:exact))Of course, creating a representation in floating points with exact arithmetic works here because we have 0.5 which is 0.1 in binary but in general, is not a good idea.julia> Rational{BigInt}(1/2)\n1//2\njulia> Rational{BigInt}(1/3)\n6004799503160661//18014398509481984\njulia> Rational{BigInt}(1/5)\n3602879701896397//18014398509481984"
},

{
    "location": "polyhedron/#Retrieving-a-representation-1",
    "page": "Polyhedron",
    "title": "Retrieving a representation",
    "category": "section",
    "text": "One can retrieve an H-representation (resp. V-representation) from a polyhedron using hrep (resp. vrep). The concrete subtype of HRepresentation (resp. VRepresentation) returned is not necessarily the same that the one used to create the polyhedron. As a rule of thumb, it is the representation the closest to the internal representation used by the library.julia> hr = hrep(poly)\njulia> typeof(hr)\nPolyhedra.LiftedHRepresentation{2,Rational{BigInt}}\njulia> hr = MixedMatHRep(hr)\njulia> typeof(hr)\nPolyhedra.MixedMatHRep{2,Rational{BigInt}}\njulia> hr.A\n3x2 Array{Rational{BigInt},2}:\n  1//1   1//1\n  1//1  -1//1\n -1//1   0//1\njulia> hr.b\n3-element Array{Rational{BigInt},1}:\n 1//1\n 0//1\n 0//1\njulia> vr = vrep(poly)\njulia> typeof(vr)\nPolyhedra.LiftedVRepresentation{2,Rational{BigInt}}\njulia> vr = MixedMatVRep(vrep)\njulia> typeof(vr)\nPolyhedra.MixedMatVRep{2,Rational{BigInt}}\njulia> vr.V\n3x2 Array{Rational{BigInt},2}:\n 1//2  1//2\n 0//1  1//1\n 0//1  0//1\n\njulia> vr.R\n0x2 Array{Rational{BigInt},2}"
},

{
    "location": "polyhedron/#Polyhedra.hrepiscomputed",
    "page": "Polyhedron",
    "title": "Polyhedra.hrepiscomputed",
    "category": "function",
    "text": "hrepiscomputed(p::Polyhedron)\n\nReturns whether the H-representation of this polyhedron has been computed.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.vrepiscomputed",
    "page": "Polyhedron",
    "title": "Polyhedra.vrepiscomputed",
    "category": "function",
    "text": "vrepiscomputed(p::Polyhedron)\n\nReturns whether the V-representation of this polyhedron has been computed.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Checking-if-a-representation-has-been-computed-1",
    "page": "Polyhedron",
    "title": "Checking if a representation has been computed",
    "category": "section",
    "text": "hrepiscomputed\nvrepiscomputed"
},

{
    "location": "polyhedron/#Polyhedra.Index",
    "page": "Polyhedron",
    "title": "Polyhedra.Index",
    "category": "type",
    "text": "Index{T, ElemT}\n\nIndex of an element of type ElemT in a Rep{T}.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.Indices",
    "page": "Polyhedron",
    "title": "Polyhedra.Indices",
    "category": "type",
    "text": "Indices{T, ElemT, RepT<:Rep{T}}\n\nIterator over the indices of the elements of type ElemT of the field rep.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidenthalfspaces",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthalfspaces",
    "category": "function",
    "text": "incidenthalfspaces(p::Polyhedron, idx)\n\nReturns the list of halfspaces incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidenthalfspaceindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthalfspaceindices",
    "category": "function",
    "text": "incidenthalfspaceindices(p::Polyhedron, idx)\n\nReturns the list of the indices of halfspaces incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentpoints",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentpoints",
    "category": "function",
    "text": "incidentpoints(p::Polyhedron, idx)\n\nReturns the list of points incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentpointindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentpointindices",
    "category": "function",
    "text": "incidentpointindices(p::Polyhedron, idx)\n\nReturns the list of the indices of points incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentrays",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentrays",
    "category": "function",
    "text": "incidentrays(p::Polyhedron, idx)\n\nReturns the list of rays incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentrayindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentrayindices",
    "category": "function",
    "text": "incidentrayindices(p::Polyhedron, idx)\n\nReturns the list of the indices of rays incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidenthyperplanes",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthyperplanes",
    "category": "function",
    "text": "incidenthyperplanes(p::Polyhedron, idx)\n\nReturns the list of hyperplanes incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidenthyperplaneindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthyperplaneindices",
    "category": "function",
    "text": "incidenthyperplaneindices(p::Polyhedron, idx)\n\nReturns the list of the indices of hyperplanes incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentlines",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentlines",
    "category": "function",
    "text": "incidentlines(p::Polyhedron, idx)\n\nReturns the list of lines incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.incidentlineindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentlineindices",
    "category": "function",
    "text": "incidentlineindices(p::Polyhedron, idx)\n\nReturns the list of the indices of lines incident to idx for the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Incidence-1",
    "page": "Polyhedron",
    "title": "Incidence",
    "category": "section",
    "text": "Elements can be accessed in a representation or polyhedron using indices and Base.get:Polyhedra.Index\nPolyhedra.IndicesThe list of indices can be obtained using, e.g., eachindex(points(rep)). For instance, the following prints all points using indicesfor pi in eachindex(points(rep))\n    @show get(rep, pi)\nendA point p (resp. ray r) is incident to an halfspace langle a x rangle le beta if langle a p rangle = beta (resp. langle a r rangle = beta).incidenthalfspaces\nincidenthalfspaceindices\nincidentpoints\nincidentpointindices\nincidentrays\nincidentrayindicesIn a polyhedron, all points and rays are incident to all hyperplanes and all halfspaces are incident to all lines. The following methods are therefore redundant, e.g. incidenthyperplanes(p, idx) is equivalent to hyperplanes(p) and incidenthyperplaneindices(p, idx) is equivalent to eachindex(hyperplanes(p)). The methods are hence only defined for consistency.incidenthyperplanes\nincidenthyperplaneindices\nincidentlines\nincidentlineindices"
},

{
    "location": "polyhedron/#Polyhedra.default_library",
    "page": "Polyhedron",
    "title": "Polyhedra.default_library",
    "category": "function",
    "text": "default_library(d::FullDim, ::Type{T}) where {T}\n\nReturns the default polyhedral library for d-dimensional polyhedron of coefficient type T.\n\nExamples\n\nTo obtain the default library for 2-dimensional polyhedra of eltype Float64, do default_library(2, Float64).\n\nGiven an StaticArrays.SVector v, to obtain a default library for points of the type of v in a type stable way, do default_library(Polyhedra.FullDim(v), eltype(v)).\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.similar_library",
    "page": "Polyhedron",
    "title": "Polyhedra.similar_library",
    "category": "function",
    "text": "similar_library(lib::Library, d::FullDim, T::Type)\n\nReturns a library that supports polyhedra of full dimension T with coefficient type T. If lib does not support it, this commonly calls default_library(d, T).\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.library",
    "page": "Polyhedron",
    "title": "Polyhedra.library",
    "category": "function",
    "text": "library(p::Polyhedron)\n\nReturns the library used by p.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.default_type",
    "page": "Polyhedron",
    "title": "Polyhedra.default_type",
    "category": "function",
    "text": "default_type(d::FullDim, ::Type{T}) where {T}\n\nReturns the default polyhedron type for d-dimensional polyhedron of coefficient type T.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.DefaultLibrary",
    "page": "Polyhedron",
    "title": "Polyhedra.DefaultLibrary",
    "category": "type",
    "text": "DefaultLibrary{T}\n\nDefault library for polyhedra of dimension larger than 1 (IntervalLibrary is the default for polyhedra of dimension 1). The library implements the bare minimum and uses the fallback implementation for all operations.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Polyhedra.IntervalLibrary",
    "page": "Polyhedron",
    "title": "Polyhedra.IntervalLibrary",
    "category": "type",
    "text": "IntervalLibrary{T}\n\nDefault library for polyhedra of dimension 1. Many aspect of polyhedral computation become trivial in one dimension. This library exploits this fact. The library is also used as a fallback for libraries that do not support 1-dimensional polyhedra (e.g. qhull). That is projecting a polyhedron using such library produces a polyhedron using IntervalLibrary.\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Base.similar",
    "page": "Polyhedron",
    "title": "Base.similar",
    "category": "function",
    "text": "similar(p::Tuple{Vararg{Polyhedra.Rep}}, d::Polyhedra.FullDim, ::Type{T}, it::Polyhedra.It{T}...)\n\nCreates a representation with a type similar to p of a polyhedron of full dimension d, element type T and initialize it with the iterators it. The type of the result will be chosen closer to the type of p[1].\n\n\n\n\n\n"
},

{
    "location": "polyhedron/#Default-libraries-1",
    "page": "Polyhedron",
    "title": "Default libraries",
    "category": "section",
    "text": "The following functions allows to select a default library:default_library\nsimilar_library\nlibrary\ndefault_typeThe following libraries serves as fallback:DefaultLibrary\nIntervalLibraryThe type and library of the polyhedron obtained after applying an operation of several polyhedra (of possibly different type and/or library) is determined by the similar function.similar"
},

{
    "location": "plot/#",
    "page": "Plot",
    "title": "Plot",
    "category": "page",
    "text": ""
},

{
    "location": "plot/#Plot-1",
    "page": "Plot",
    "title": "Plot",
    "category": "section",
    "text": "Polyhedra contains utilities to visualize either a 2-dimensional or a 3-dimensional polyhedron, see Polyhedron for how to construct a polyhedron, e.g. from its H- or V-representation."
},

{
    "location": "plot/#D-plotting-with-Plots-1",
    "page": "Plot",
    "title": "2D plotting with Plots",
    "category": "section",
    "text": "A 2-dimensional polytope, i.e. bounded polyhedron, can be visualized with Plots. Suppose for instance that we want to visualize the polyhedron having the following H-representation:julia> using Polyhedra\n\njulia> h = HalfSpace([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0)\nH-representation Polyhedra.Intersection{Int64,Array{Int64,1},Int64}:\n3-element iterator of HalfSpace{Int64,Array{Int64,1}}:\n HalfSpace([1, 1], 1)\n HalfSpace([-1, 0], 0)\n HalfSpace([0, -1], 0)The H-representation cannot be given to Plots directly, it first need to be transformed into a polyhedron:julia> p = polyhedron(h)\nPolyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},Int64},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},Int64}}:\n3-element iterator of HalfSpace{Rational{BigInt},Array{Rational{BigInt},1}}:\n HalfSpace(Rational{BigInt}[1//1, 1//1], 1//1)\n HalfSpace(Rational{BigInt}[-1//1, 0//1], 0//1)\n HalfSpace(Rational{BigInt}[0//1, -1//1], 0//1)The polyhedron can be given to Plots as followsjulia> using Plots\n\njulia> plot(p)See Polyhedral Function and 3D Plotting a projection of the 4D permutahedron for example notebooks."
},

{
    "location": "plot/#Polyhedra.Mesh",
    "page": "Plot",
    "title": "Polyhedra.Mesh",
    "category": "type",
    "text": "struct Mesh{N, T, PT <: Polyhedron{T}} <: GeometryTypes.GeometryPrimitive{N, T}\n    polyhedron::PT\nend\n\nMesh wrapper type that inherits from GeometryPrimitive to be used for plotting a polyhedron. Note that Mesh(p) is type unstable but one can use Mesh{3}(p) instead if it is known that p is defined in a 3-dimensional space.\n\n\n\n\n\n"
},

{
    "location": "plot/#D-plotting-with-Plots-2",
    "page": "Plot",
    "title": "3D plotting with Plots",
    "category": "section",
    "text": "A 3-dimensional polyhedron can be visualized with either MeshCat or Makie. Unbounded polyhedron are supported by truncating the polyhedron into a polytope and not triangularizing the faces in the directions of unbounded rays.Suppose for instance that we want to visualize the polyhedron having the following H-representation:julia> using Polyhedra\n\njulia> v = convexhull([0, 0, 0]) + conichull([1, 0, 0], [0, 1, 0], [0, 0, 1])\nV-representation Polyhedra.Hull{Int64,Array{Int64,1},Int64}:\n1-element iterator of Array{Int64,1}:\n [0, 0, 0],\n3-element iterator of Ray{Int64,Array{Int64,1}}:\n Ray([1, 0, 0])\n Ray([0, 1, 0])\n Ray([0, 0, 1])The V-representation cannot be given to MeshCat or Makie directly, it first need to be transformed into a polyhedron:julia> p = polyhedron(v)\nPolyhedron DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},Int64},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},Int64}}:\n1-element iterator of Array{Rational{BigInt},1}:\n Rational{BigInt}[0//1, 0//1, 0//1],\n3-element iterator of Ray{Rational{BigInt},Array{Rational{BigInt},1}}:\n Ray(Rational{BigInt}[1//1, 0//1, 0//1])\n Ray(Rational{BigInt}[0//1, 1//1, 0//1])\n Ray(Rational{BigInt}[0//1, 0//1, 1//1])Then, we need to create a mess from the polyhedron:julia> m = Polyhedra.Mesh(p)\nPolyhedra.Mesh{3,Rational{BigInt},DefaultPolyhedron{Rational{BigInt},Polyhedra.Intersection{Rational{BigInt},Array{Rational{BigInt},1},Int64},Polyhedra.Hull{Rational{BigInt},Array{Rational{BigInt},1},Int64}}}(convexhull([0//1, 0//1, 0//1]) + convexhull(Ray(Rational{BigInt}[1//1, 0//1, 0//1]), Ray(Rational{BigInt}[0//1, 1//1, 0//1]), Ray(Rational{BigInt}[0//1, 0//1, 1//1])))Polyhedra.MeshThe polyhedron can be plotted with MeshCat as followsjulia> using MeshCat\n\njulia> vis = Visualizer()\n\njulia> setobject!(vis, m)\n\njulia> open(vis)To plot it in a notebook, replace open(vis) with IJuliaCell(vis).To plot it with Makie instead, you can use for instance mesh or wireframe.julia> import Makie\n\njulia> Makie.mesh(m, color=:blue)\n\njulia> Makie.wireframe(m)See 3D Plotting a projection of the 4D permutahedron for an example notebook."
},

{
    "location": "redundancy/#",
    "page": "Containment/Redundancy",
    "title": "Containment/Redundancy",
    "category": "page",
    "text": ""
},

{
    "location": "redundancy/#Containment/Redundancy-1",
    "page": "Containment/Redundancy",
    "title": "Containment/Redundancy",
    "category": "section",
    "text": ""
},

{
    "location": "redundancy/#Base.in",
    "page": "Containment/Redundancy",
    "title": "Base.in",
    "category": "function",
    "text": "in(p::VRepElement, h::HRepElement)\n\nReturns whether p is in h. If h is an hyperplane, it returns whether langle a x rangle approx beta. If h is an halfspace, it returns whether langle a x rangle le beta.\n\nin(p::VRepElement, h::HRep)\n\nReturns whether p is in h, e.g. in all the hyperplanes and halfspaces supporting h.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Base.issubset",
    "page": "Containment/Redundancy",
    "title": "Base.issubset",
    "category": "function",
    "text": "issubset(p::Rep, h::HRepElement)\n\nReturns whether p is a subset of h, i.e. whether h supports the polyhedron p.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.ininterior",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.ininterior",
    "category": "function",
    "text": "ininterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the interior of h. If h is an hyperplane, it always returns false. If h is an halfspace langle a x rangle leq beta, it returns whether p is in the open halfspace langle a x rangle  beta\n\nininterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the interior of h, e.g. in the interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.inrelativeinterior",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.inrelativeinterior",
    "category": "function",
    "text": "inrelativeinterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the relative interior of h. If h is an hyperplane, it is equivalent to p in h since the relative interior of an hyperplane is itself. If h is an halfspace, it is equivalent to ininterior(p, h).\n\ninrelativeinterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the relative interior of h, e.g. in the relative interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.support_function",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.support_function",
    "category": "function",
    "text": "supportfunction(h::AbstractVector, rep::Rep, solver=Polyhedra.linearobjective_solver(p))\n\nReturn the value of the support function of rep at h. See [Section 13, R15] for more details.\n\n[R15] Rockafellar, R.T. Convex analysis. Princeton university press, 2015.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Containment-1",
    "page": "Containment/Redundancy",
    "title": "Containment",
    "category": "section",
    "text": "in\nissubset\nininterior\ninrelativeinterior\nsupport_function"
},

{
    "location": "redundancy/#Polyhedra.detecthlinearity!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.detecthlinearity!",
    "category": "function",
    "text": "detecthlinearity!(p::VRep)\n\nDetects all the hyperplanes contained in the H-representation and remove all redundant hyperplanes.\n\nExamples\n\nThe representation\n\nh = HalfSpace([1, 1], 1]) ∩ HalfSpace([-1, -1], -1)\n\ncontains the hyperplane HyperPlane([1, 1], 1).\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.detectvlinearity!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.detectvlinearity!",
    "category": "function",
    "text": "detectvlinearity!(p::VRep)\n\nDetects all the lines contained in the V-representation and remove all redundant lines.\n\nExamples\n\nThe representation\n\nv = conichull([1, 1], [-1, -1])\n\ncontains the line Line([1, 1]).\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.dim",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.dim",
    "category": "function",
    "text": "dim(h::HRep, current=false)\n\nReturns the dimension of the affine hull of the polyhedron. That is the number of non-redundant hyperplanes that define it. If current is true then it simply returns the dimension according the current number of hyperplanes, assuming that the H-linearity has already been detected. Otherwise, it first calls detecthlinearity!.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Linearity-1",
    "page": "Containment/Redundancy",
    "title": "Linearity",
    "category": "section",
    "text": "detecthlinearity!\ndetectvlinearity!\ndim"
},

{
    "location": "redundancy/#Polyhedra.removeduplicates",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removeduplicates",
    "category": "function",
    "text": "removeduplicates(rep::Representation)\n\nRemoves the duplicates in the Representation.\n\nIn an H-representation, it removes the redundant hyperplanes and it remove an halfspace when it is equal to another halfspace in the affine hull. For instance, HalfSpace([1, 1], 1) is equal to HalfSpace([1, 0], 0) in the affine hull generated by HyperPlane([0, 1], 1]).\nIn a V-representation, it removes the redundant lines and it remove a point (resp. ray) when it is equal to another point (resp. ray) in the line hull. For instance, in the line hull generated by Line([0, 1]), [1, 1] is equal to [1, 0] and Ray([2, 2]) is equal to Ray([1, 0]).\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Duplicates-1",
    "page": "Containment/Redundancy",
    "title": "Duplicates",
    "category": "section",
    "text": "removeduplicates"
},

{
    "location": "redundancy/#Polyhedra.isredundant",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.isredundant",
    "category": "function",
    "text": "isredundant(p::Rep, idx::Index; strongly=false)\n\nReturn a Bool indicating whether the element with index idx can be removed without changing the polyhedron represented by p. If strongly is true,\n\nif idx is an H-representation element h, it returns true only if no V-representation element of p is in the hyperplane of h.\nif idx is a V-representation element v, it returns true only if v is in the relative interior of p.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.removehredundancy!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removehredundancy!",
    "category": "function",
    "text": "removehredundancy!(p::HRep)\n\nRemoves the elements of the H-representation of p that can be removed without changing the polyhedron represented by p. That is, it only keeps the halfspaces corresponding to facets of the polyhedron.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.removevredundancy",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removevredundancy",
    "category": "function",
    "text": "removevredundancy(vr::VRepresentation)\n\nReturn a V-representation of the polyhedron represented by vr all the elements of vr except the redundant ones, i.e. the elements that can be expressed as convex combination of other ones.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Polyhedra.removevredundancy!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removevredundancy!",
    "category": "function",
    "text": "removevredundancy!(p::VRep; strongly=false)\n\nRemoves the elements of the V-representation of p that can be removed without changing the polyhedron represented by p. That is, it only keeps the extreme points and rays. This operation is often called \"convex hull\" as the remaining points are the extreme points of the convex hull of the initial set of points. If strongly=true, weakly redundant points, i.e., points that are not extreme but are not in the relative interior either, may be kept.\n\n\n\n\n\n"
},

{
    "location": "redundancy/#Redundancy-1",
    "page": "Containment/Redundancy",
    "title": "Redundancy",
    "category": "section",
    "text": "isredundant\nremovehredundancy!\nremovevredundancy\nremovevredundancy!"
},

{
    "location": "projection/#",
    "page": "Projection/Elimination",
    "title": "Projection/Elimination",
    "category": "page",
    "text": ""
},

{
    "location": "projection/#Polyhedra.FourierMotzkin",
    "page": "Projection/Elimination",
    "title": "Polyhedra.FourierMotzkin",
    "category": "type",
    "text": "FourierMotzkin\n\nComputation of the projection by computing the H-representation and applying the Fourier-Motzkin elimination algorithm to it.\n\n\n\n\n\n"
},

{
    "location": "projection/#Polyhedra.BlockElimination",
    "page": "Projection/Elimination",
    "title": "Polyhedra.BlockElimination",
    "category": "type",
    "text": "BlockElimination\n\nComputation of the projection by computing the H-representation and applying the block elimination algorithm to it.\n\n\n\n\n\n"
},

{
    "location": "projection/#Polyhedra.ProjectGenerators",
    "page": "Projection/Elimination",
    "title": "Polyhedra.ProjectGenerators",
    "category": "type",
    "text": "ProjectGenerators\n\nComputation of the projection by computing the V-representation and projecting them.\n\n\n\n\n\n"
},

{
    "location": "projection/#Polyhedra.eliminate",
    "page": "Projection/Elimination",
    "title": "Polyhedra.eliminate",
    "category": "function",
    "text": "eliminate(p::Polyhedron, delset, algo::EliminationAlgorithm)\n\nEliminate the dimensions in delset by projecting the polyhedron onto the remaining dimension.\n\n\n\n\n\n"
},

{
    "location": "projection/#Polyhedra.project",
    "page": "Projection/Elimination",
    "title": "Polyhedra.project",
    "category": "function",
    "text": "project(p::Polyhedron, pset, algo)\n\nEquivalent to eliminate(p, setdiff(1:fulldim(p), pset), algo).\n\n\n\n\n\n"
},

{
    "location": "projection/#Polyhedra.fixandeliminate",
    "page": "Projection/Elimination",
    "title": "Polyhedra.fixandeliminate",
    "category": "function",
    "text": "fixandeliminate(p::HRep{T}, I, v)\n\nFix the variables with indices in I to the corresponding value in v. This is equivalent to doing the following:\n\nfunction ei(i)\n    a = zeros(T, fulldim(p))\n    a[i] = one(T)\n    a\nend\neliminate(p ∩ HyperPlane(ei(I[1]), v[1]) ∩ ... ∩ HyperPlane(ei(I[n]), v[n]))\n\nwhere n is the length of I (and v), but it is much more efficient. The code above does a polyhedral projection while this function simply replaces each halfspace ⟨a, x⟩ ≤ β (resp. each hyperplane ⟨a, x⟩ = β) by the halfspace ⟨a_J, x⟩ ≤ β - ⟨a_I, v⟩ (resp. the hyperplane ⟨a_J, x⟩ = β - ⟨a_I, v⟩) where J = setdiff(1:fulldim(p), I).\n\n\n\n\n\n"
},

{
    "location": "projection/#Projection/Elimination-1",
    "page": "Projection/Elimination",
    "title": "Projection/Elimination",
    "category": "section",
    "text": "Consider the polyhedron created in the beginning of this section. As a reminder, it represents the following H-representation:beginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*One can verify that for any 0 leq x_2 leq 1, there exists a value x_1 such that (x_1 x_2) is in this polyhedron. This means that the H-representation obtained by eliminating x_1 is:beginalign*\n  x_1  leq 1 \n  x_1  geq 0\nendalign*where x_1 in the H-representation above represents x_2 in the previous one. This can be obtained as followsjulia> poly_x2 = eliminate(poly, [1])\njulia> hrep(poly_x2)\nH-representation\nbegin\n 2 2 rational\n 1//1 -1//1\n 0//1 1//1\nendThere is two methods of computing the elimination implemented in CDDLib: Fourier-Motzkin elimination and block elimination. As written by K. Fukuda in CDD\'s documentation, \"[Block elimination] might be a faster way to eliminate variables than the repeated [Fourier-Motzkin elimination] when the number of variables to eliminate is large\". You can specify the method to use as a third argument, e.g. eliminate(poly, [1], FourierMotzkin()), eliminate(poly, [1], BlockElimination()). A third method can be chosen: ProjectGenerators. It computes the V-representation and then project each of its elements. This is the method of choice when the V-representation is already computed.FourierMotzkin\nBlockElimination\nProjectGeneratorsIf nothing is specified as in the block of code above, the behavior depends on the polyhedral library. If neither Fourier-Motzkin nor block elimination is implemented or if the V-representation is already computed then :ProjectGenerators is chosen. Otherwise, Polyhedra lets the library decide. In CDDLib, :FourierMotzkin is chosen when only the last dimension needs to be eliminated and :BlockElimination is chosen otherwise. Note that CDDLib only supports projecting the last trailing dimensions.eliminate\nproject\nfixandeliminate"
},

{
    "location": "optimization/#",
    "page": "Optimization",
    "title": "Optimization",
    "category": "page",
    "text": ""
},

{
    "location": "optimization/#Base.isempty",
    "page": "Optimization",
    "title": "Base.isempty",
    "category": "function",
    "text": "isempty(p::Rep, solver::JuMP.OptimizerFactory=Polyhedra.linear_objective_solver(p))\n\nCheck whether the polyhedron p is empty by using the solver solver.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Polyhedra.VRepOptimizer",
    "page": "Optimization",
    "title": "Polyhedra.VRepOptimizer",
    "category": "type",
    "text": "VRepOptimizer{T} <: AbstractPolyhedraOptimizer{T}\n\nLinear Programming solver using the V-representation of the feasible set to find the optimal solution.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Polyhedra.default_solver",
    "page": "Optimization",
    "title": "Polyhedra.default_solver",
    "category": "function",
    "text": "default_solver(p::Rep)\n\nReturns a default linear programming solver for the polyhedron p (e.g. CDD has an internal solver which is used by default).\n\n\n\n\n\n"
},

{
    "location": "optimization/#Polyhedra.linear_objective_solver",
    "page": "Optimization",
    "title": "Polyhedra.linear_objective_solver",
    "category": "function",
    "text": "linear_objective_solver(p::Rep, solver::Union{Nothing, JuMP.OptimizerFactory}=default_solver(p))\n\nReturn the solver to use for optimizing a linear objective over the polyhedron p, i.e.\n\nmodel = Model(solver)\nx = @variable(model, [1:fulldim(p)])\n@constraint(model, x in p)\n@objective(model, c ⋅ x)\n\nfor some vector c.\n\nBy default, if the V-representation of p has been computed, it returns VRepOptimizer(), otherwise, it returns solver.\n\nIf the problem has constraints different to x in p, use default_solver(p) instead as the fact that the V-representation of p has been computed does not help.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Optimization-1",
    "page": "Optimization",
    "title": "Optimization",
    "category": "section",
    "text": "A polyhedron can represents the feasible set of an optimization program. The program is infeasible when the polyhedron is empty.isemptyIf the V-representation of the polyhedron has been computed, it can be used to solve the linear program.VRepOptimizerOtherwise, any programming solver implementing the MathOptInterface interface can be used. See here for a list of available solvers.Polyhedra.default_solver\nPolyhedra.linear_objective_solver"
},

{
    "location": "optimization/#Polyhedra.PolyhedraToLPBridge",
    "page": "Optimization",
    "title": "Polyhedra.PolyhedraToLPBridge",
    "category": "type",
    "text": "PolyhedraToLPBridge{T}\n\nThe PolyhedraToLPBridge converts a constraint VF-in-PolyhedraOptSet into the constraints F-in-EqualTo for the hyperplanes and F-to-LessThan for halfspaces.\n\n\n\n\n\n"
},

{
    "location": "optimization/#Using-a-polyhedron-for-in-an-optimization-model-1",
    "page": "Optimization",
    "title": "Using a polyhedron for in an optimization model",
    "category": "section",
    "text": "A polyhedron or representation can be used in the constraint of a JuMP model. For instance, consider the 1-simplex:julia> using Polyhedra\n\njulia> simplex = HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0) ∩ HyperPlane([1, 1], 1)\nH-representation Polyhedra.Intersection{Int64,Array{Int64,1},Int64}:\n1-element iterator of HyperPlane{Int64,Array{Int64,1}}:\n HyperPlane([1, 1], 1),\n2-element iterator of HalfSpace{Int64,Array{Int64,1}}:\n HalfSpace([-1, 0], 0)\n HalfSpace([0, -1], 0)and the following JuMP model with two variablesjulia> using JuMP\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, λ[1:2])\n2-element Array{VariableRef,1}:\n λ[1]\n λ[2]The variables can be constrained to belong to the simplex as follows:julia> @constraint(model, λ in simplex)\n[λ[1], λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))but a vector of affine or quadratic expressions can also be constrained to belong to the simplex:julia> A = [1  1\n            1 -1]\n2×2 Array{Int64,2}:\n 1   1\n 1  -1\n\njulia> @constraint(model, A * λ in simplex)\n[λ[1] + λ[2], λ[1] - λ[2]] ∈ Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}(HyperPlane([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0))We can verify that the model contains both constraints:julia> model\nA JuMP Model\nFeasibility problem with:\nVariables: 2\n`Array{JuMP.VariableRef,1}`-in-`Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}`: 1 constraint\n`Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},1}`-in-`Polyhedra.PolyhedraOptSet{Int64,Polyhedra.Intersection{Int64,Array{Int64,1},Int64}}`: 1 constraint\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\nNames registered in the model: λWhen the model is solved, the constraint is automatically transformed into appropriate constraints if the optimizer does not support consraints with the set Polyhedra.PolyhedraOptSet.julia> import GLPK\n\njulia> optimize!(model, with_optimizer(GLPK.Optimizer))\n\njulia> termination_status(model)\nOPTIMAL::TerminationStatusCode = 1\n\njulia> value.(λ)\n2-element Array{Float64,1}:\n 0.5\n 0.5For instance, GLPK, does not support Polyhedra.PolyhedraOptSet constraints but supports MOI.EqualTo{Float64} and MOI.LessThan{Float64}. The polyhedral constraints are therefore bridged into several MOI.EqualTo{Float64} and MOI.LessThan{Float64} constraints using the following constraint bridge:Polyhedra.PolyhedraToLPBridgeSee Polyhedral Function for an example notebook."
},

{
    "location": "optimization/#Creating-a-polyhedron-from-the-feasible-set-of-a-JuMP-model-1",
    "page": "Optimization",
    "title": "Creating a polyhedron from the feasible set of a JuMP model",
    "category": "section",
    "text": "A typical application of polyhedral computation is the computation of the set of extreme points and rays of the feasible set of an optimization problem. This comes from the fact that given a minimization of a concave function (or maximization of a convex function) on a convex feasible set (e.g. Linear Programming), we are either in the following three situations:The feasible set is empty, i.e. the problem is infeasible.\nAn extreme ray is optimal, i.e. the problem is unbounded (or it may also be bounded if the objective is constant along the ray).\nAn extreme point is optimal.A JuMP model is treated by polyhedron just like any H-representation. For example, the hypercube of dimension n can be created as followsm = Model()\n@variable(m, 0 ≤ x[1:n] ≤ 1)\n\npoly = polyhedron(m, CDDLib.Library(:exact))"
},

{
    "location": "utilities/#",
    "page": "Utilities",
    "title": "Utilities",
    "category": "page",
    "text": ""
},

{
    "location": "utilities/#Utilities-1",
    "page": "Utilities",
    "title": "Utilities",
    "category": "section",
    "text": ""
},

{
    "location": "utilities/#Base.:*",
    "page": "Utilities",
    "title": "Base.:*",
    "category": "function",
    "text": "*(p1::Rep, p2::Rep)\n\nCartesian product between the polyhedra p1 and p2.\n\n\n\n\n\n*(P::Union{AbstractMatrix, UniformScaling}, p::VRep)\n\nTransform the polyhedron represented by p into P p by transforming each element of the V-representation (points, symmetric points, rays and lines) x into P x.\n\n\n\n\n\n*(α::Number, p::Rep)\n\nTransform the polyhedron represented by p into alpha p by transforming each element of the V-representation (points, symmetric points, rays and lines) x into alpha x.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Base.:\\",
    "page": "Utilities",
    "title": "Base.:\\",
    "category": "function",
    "text": "(P::Union{AbstractMatrix, UniformScaling}, p::HRep)\n\nTransform the polyhedron represented by p into P^-1 p by transforming each halfspace langle a x rangle le beta into langle P^top a x rangle le beta and each hyperplane langle a x rangle = beta into langle P^top a x rangle = beta.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Base.:/",
    "page": "Utilities",
    "title": "Base.:/",
    "category": "function",
    "text": "/(p::HRep, P::Union{AbstractMatrix, UniformScaling})\n\nTransform the polyhedron represented by p into P^-T p by transforming each halfspace langle a x rangle le beta into langle P a x rangle le beta and each hyperplane langle a x rangle = beta into langle P a x rangle = beta.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Base.intersect",
    "page": "Utilities",
    "title": "Base.intersect",
    "category": "function",
    "text": "intersect(P1::HRep, P2::HRep)\n\nTakes the intersection of P1 and P2  x  x in P_1 x in P_2 . It is very efficient between two H-representations or between two polyhedron for which the H-representation has already been computed. However, if P1 (resp. P2) is a polyhedron for which the H-representation has not been computed yet, it will trigger a representation conversion which is costly. See the Polyhedral Computation FAQ for a discussion on this operation.\n\nThe type of the result will be chosen closer to the type of P1. For instance, if P1 is a polyhedron (resp. H-representation) and P2 is a H-representation (resp. polyhedron), intersect(P1, P2) will be a polyhedron (resp. H-representation). If P1 and P2 are both polyhedra (resp. H-representation), the resulting polyhedron type (resp. H-representation type) will be computed according to the type of P1. The coefficient type however, will be promoted as required taking both the coefficient type of P1 and P2 into account.\n\n\n\n\n\nintersect(v::VRepresentation{T}, h::HRepElement)\n\nCompute the intersection of v with an halfspace or hyperplane h. The method used by default is to keep the V-representation element of v that are in h and add new ones generated as the intersection between the hyperplane defining h and the segment between two adjacent V-representation elements of v that are in either sides of the hyperplane. See Lemma 3 of [FP96] for more detail on the method.\n\n[FP96] Fukuda, K. and Prodon, A. Double description method revisited Combinatorics and computer science, Springer, 1996, 91-111\n\n\n\n\n\n"
},

{
    "location": "utilities/#Base.intersect!",
    "page": "Utilities",
    "title": "Base.intersect!",
    "category": "function",
    "text": "intersect!(p::HRep, h::Union{HRepresentation, HRepElement})\n\nSame as intersect except that p is modified to be equal to the intersection.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.convexhull",
    "page": "Utilities",
    "title": "Polyhedra.convexhull",
    "category": "function",
    "text": "convexhull(P1::VRep, P2::VRep)\n\nTakes the convex hull of P1 and P2  lambda x + (1-lambda) y  x in P_1 y in P_2 . It is very efficient between two V-representations or between two polyhedron for which the V-representation has already been computed. However, if P1 (resp. P2) is a polyhedron for which the V-representation has not been computed yet, it will trigger a representation conversion which is costly.\n\nThe type of the result will be chosen closer to the type of P1. For instance, if P1 is a polyhedron (resp. V-representation) and P2 is a V-representation (resp. polyhedron), convexhull(P1, P2) will be a polyhedron (resp. V-representation). If P1 and P2 are both polyhedra (resp. V-representation), the resulting polyhedron type (resp. V-representation type) will be computed according to the type of P1. The coefficient type however, will be promoted as required taking both the coefficient type of P1 and P2 into account.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.convexhull!",
    "page": "Utilities",
    "title": "Polyhedra.convexhull!",
    "category": "function",
    "text": "convexhull!(p1::VRep, p2::VRep)\n\nSame as convexhull except that p1 is modified to be equal to the convex hull.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.translate",
    "page": "Utilities",
    "title": "Polyhedra.translate",
    "category": "function",
    "text": "translate(p::Polyhedra.Rep, v::AbstractVector)\n\nComputes translation of the polyhedron p with the vector v. That is, computes\n\n x + v mid x in p \n\nBy default, if the H-representation, it simply translates every hyperplanes and halfspace, otherwise, it translates every points of the V-representation. That is, this operation can be achieved both in the H-representation and V-representation hence does not trigger any representation conversion.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Operations-1",
    "page": "Utilities",
    "title": "Operations",
    "category": "section",
    "text": "*\n\\\n/\nintersect\nintersect!\nconvexhull\nconvexhull!\ntranslate"
},

{
    "location": "utilities/#Polyhedra.volume",
    "page": "Utilities",
    "title": "Polyhedra.volume",
    "category": "function",
    "text": "volume(p::Polyhedron{T}) where {T}\n\nReturns the fulldim(p)-dimensional hyper-volume of the polyhedron p. Returns Inf or -one(T) if it is infinite depending on whether the type T has an infinite value.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.surface",
    "page": "Utilities",
    "title": "Polyhedra.surface",
    "category": "function",
    "text": "surface(p::Polyhedron{T}) where {T}\n\nReturns the fulldim(p)-1-dimensional hyper-volume of the surface of the polyhedron p. Returns Inf or -one(T) if it is infinite depending on whether the type T has an infinite value.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Volume-1",
    "page": "Utilities",
    "title": "Volume",
    "category": "section",
    "text": "volume\nsurface"
},

{
    "location": "utilities/#Polyhedra.chebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.chebyshevcenter",
    "category": "function",
    "text": "chebyshevcenter(p::Rep[, solver])\n\nIf p is a H-representation or is a polyhedron for which the H-representation has already been computed, calls hchebyshevcenter, otherwise, call vchebyshevcenter.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.hchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.hchebyshevcenter",
    "category": "function",
    "text": "hchebyshevcenter(p::HRep[, solver])\n\nReturn a tuple with the center and radius of the largest euclidean ball contained in the polyhedron p. Throws an error if the polyhedron is empty or if the radius is infinite.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.vchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.vchebyshevcenter",
    "category": "function",
    "text": "vchebyshevcenter(p::VRep[, solver])\n\nReturn a tuple with the center and radius of the smallest euclidean ball containing the polyhedron p. Throws an error if the polyhedron is empty or if the radius is infinite (i.e. p is not a polytope, it contains rays).\n\n\n\n\n\n"
},

{
    "location": "utilities/#Chebyshev-center-1",
    "page": "Utilities",
    "title": "Chebyshev center",
    "category": "section",
    "text": "chebyshevcenter\nhchebyshevcenter\nvchebyshevcenter"
},

{
    "location": "utilities/#Polyhedra.@subrepelem",
    "page": "Utilities",
    "title": "Polyhedra.@subrepelem",
    "category": "macro",
    "text": "The representation rep contain the elements elem inside a representation in the field field.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.@norepelem",
    "page": "Utilities",
    "title": "Polyhedra.@norepelem",
    "category": "macro",
    "text": "The representation rep does not contain any elem.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Polyhedra.@vecrepelem",
    "page": "Utilities",
    "title": "Polyhedra.@vecrepelem",
    "category": "macro",
    "text": "The representation rep contain the elements elem inside a vector in the field field.\n\n\n\n\n\n"
},

{
    "location": "utilities/#Defining-new-representation-1",
    "page": "Utilities",
    "title": "Defining new representation",
    "category": "section",
    "text": "The following macros make it easy to define new representations:Polyhedra.@subrepelem\nPolyhedra.@norepelem\nPolyhedra.@vecrepelem"
},

]}
