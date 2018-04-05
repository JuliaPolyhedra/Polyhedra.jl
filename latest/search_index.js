var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Polyhedra-–-Manipulation-of-Polyhedra-in-Julia-1",
    "page": "Index",
    "title": "Polyhedra –- Manipulation of Polyhedra in Julia",
    "category": "section",
    "text": "Polyhedra is a package for polyhedra manipulations in Julia. It provides an unified interface for Polyhedra Manipulation Libraries such as CDDLib.jl, LRSLib.jl and QHull.Polyhedra can either be represented by a set of linear inequalities or by vertices and rays. In the first case, the points of the polyhedron are the points that satisfies all the inequalities and in the second case they are the points that can be expressed as a convex combination of the vertices plus a conic combination of the rays. The manipulations that Polyhedra can perform includeProjection: Projection of a polyhedron on a lower dimensional space, e.g. Fourier-Motzkin elimination.\nChanging the Representation\nVertex enumeration problem: Computing the extremal vertices and rays from an inequality representation\nConvex hull problem: Computing a set of linear inequalities describing the polyhedron from a vertex/ray representation\nRemoval of redundant inequalities or redundant vertices/rays.\nPlotting of 2D polyhedra using Plots\nDecomposition of 3D polyhedra into vertices and triangular faces, enabling easy visualization of 3D polyhedra using DrakeVisualizer or GLVisualize.Depending on the library, these manipulation can either be in floating point or exact rational arithmetic.Each operation has a default fallback implementation which is used in case the library does not support it. Polyhedra also includes a default library which does not implement anything, hence using every fallback.Polyhedra remains under active development, and we welcome your feedback, suggestions, and bug reports."
},

{
    "location": "index.html#Installing-Polyhedra-1",
    "page": "Index",
    "title": "Installing Polyhedra",
    "category": "section",
    "text": "If you are familiar with Julia you can get started quickly by using the package manager to install Polyhedrajulia> Pkg.add(\"Polyhedra\")And a Polyhedra Manipulation Library, e.g.julia> Pkg.add(\"CDDLib\")"
},

{
    "location": "index.html#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\", \"representation.md\", \"polyhedron.md\", \"redundancy.md\", \"projection.md\", \"utilities.md\"]\nDepth = 2"
},

{
    "location": "installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "This section shows how to install Julia, Polyhedra and a Polyhedra Manipulation Library of your choice."
},

{
    "location": "installation.html#Getting-Julia-1",
    "page": "Installation",
    "title": "Getting Julia",
    "category": "section",
    "text": "The first step is to install Julia. Polyhedra supports Julia v0.4 to Julia v0.6 but the latest version only supports Julia v0.6 (the command Pkg.add will automatically install the latest compatible version of Polyhedra). Download links and more detailed instructions are available on the Julia website."
},

{
    "location": "installation.html#Getting-Polyhedra-1",
    "page": "Installation",
    "title": "Getting Polyhedra",
    "category": "section",
    "text": "Open a Julia console (e.g. enter julia at the command line) and writejulia> Pkg.add(\"Polyhedra\")To start using Polyhedra, you can now just writejulia> using PolyhedraPolyhedra includes a default library supporting every operation but external library can also be used. See the next section on installing a library."
},

{
    "location": "installation.html#Getting-Libraries-1",
    "page": "Installation",
    "title": "Getting Libraries",
    "category": "section",
    "text": "Many C libraries are are available for manipulating Polyhedra. Some of them works with floating point arithmetic and some of them can do the computation exactly using rational arithmetic and multiple precision libraries such as GMP. Julia also natively support Rational arithmetic using multiple precision libraries and of course floating point arithmetic. That makes the use of both arithmetic very easy and transparent.The following table provides a list of Polyhedra Manipulation Libraries. When they have a Julia library implementing the interface of Polyhedra.jl then the \"Library\" column shows the name of the library.Solver Julia Package Library License Exact Rational Floating point\ncdd CDDLib.jl CDDLibrary() GPL X X\nConvexHull ConvexHull.jl ConvexHullLib() MIT X \nlrs LRSLib.jl LRSLibrary() GPL X X\nqhull QHull.jl QHullLib()   X\nCHull2d CHull2d.jl  MIT X X\nNewPolka None  GPL X \nParma Polyhedra Library None  GPL X \npd None  GPL X \nporta None  GPL X (overflow !) Please let me know if you plan to write a new wrapper (or an implementation in pure Julia). Since libraries use different algorithms, no library is better for every problem; here and here are comparisons."
},

{
    "location": "representation.html#",
    "page": "Representation",
    "title": "Representation",
    "category": "page",
    "text": ""
},

{
    "location": "representation.html#Polyhedra.fulldim",
    "page": "Representation",
    "title": "Polyhedra.fulldim",
    "category": "function",
    "text": "fulldim(rep::Rep)\n\nReturns the dimension of the space in which the representation is defined. That is, a straight line in a 3D space has fulldim 3.\n\n\n\n"
},

{
    "location": "representation.html#MultivariatePolynomials.coefficienttype",
    "page": "Representation",
    "title": "MultivariatePolynomials.coefficienttype",
    "category": "function",
    "text": "coefficienttype(rep::Rep)\n\nReturns the type of the coefficients used in the representation of rep.\n\n\n\n"
},

{
    "location": "representation.html#Representation-1",
    "page": "Representation",
    "title": "Representation",
    "category": "section",
    "text": "Polyhedra can be described in 2 different ways.H-representation: As the intersection of finitely many halfspaces given by its facets.\nV-representation: As the convex hull of its vertices + the conic hull of its rays where \'+\' is the Minkowski sum.In Polyhedra.jl, those representations are given the respective abstract types HRepresentation and VRepresentation which are themself subtypes of Representation.These functions can be called on both H-representation and V-representationfulldim\ncoefficienttype"
},

{
    "location": "representation.html#Polyhedra.HalfSpace",
    "page": "Representation",
    "title": "Polyhedra.HalfSpace",
    "category": "type",
    "text": "struct HalfSpace{N, T, AT} <: HRepElement{N, T}\n    a::AT\n    β::T\nend\n\nAn halfspace defined by the set of points x such that langle a x rangle le beta.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.HyperPlane",
    "page": "Representation",
    "title": "Polyhedra.HyperPlane",
    "category": "type",
    "text": "struct HyperPlane{N, T, AT} <: HRepElement{N, T}\n    a::AT\n    β::T\nend\n\nAn hyperplane defined by the set of points x such that langle a x rangle = beta.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hrep",
    "page": "Representation",
    "title": "Polyhedra.hrep",
    "category": "function",
    "text": "hrep(p::Polyhedron)\n\nReturns an H-representation for the polyhedron p.\n\n\n\nhrep(hyperplanes::HyperPlaneIt)\n\nCreates an affine space from the list of hyperplanes hyperplanes.\n\nExamples\n\nhrep([HyperPlane([0, 1, 0], 1), HyperPlane([0, 0, 1], 0)])\n\ncreates the 1-dimensional affine subspace containing all the points (x_1 0 0), i.e. the x_1-axis.\n\nhrep([HyperPlane([1, 1], 1), HyperPlane([1, 0], 0)])\n\ncreates the 0-dimensional affine subspace only containing the point (0 1).\n\n\n\nhrep(hyperplanes::HyperPlaneIt, halfspaces::HalfSpaceIt)\n\nCreates an H-representation for the polyhedron equal to the intersection of the hyperplanes hyperplanes and halfspaces halfspaces.\n\nExamples\n\nFor instance, the simplex\n\nbeginalign*\n  x_1 + x_2 = 1 \n  x_1 geq 0 \n  x_2 geq 0\nendalign*\n\ncan be created as follows:\n\nhrep([HalfSpace([-1, 0], 0)], [HyperPlane([1, 1], 1), HalfSpace([0, -1], 0)])\n\n\n\nhrep(halfspaces::HalfSpaceIt)\n\nCreates an H-representation for the polyhedron equal to the intersection of the halfspaces halfspaces.\n\nExamples\n\nFor instance, the polytope\n\neginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nndalign*\n\ncan be created as follows:\n\nhrep([HalfSpace([1, 1], 1), HalfSpace([1, -1], 0), HalfSpace([-1, 0], 0)])\n\n\n\nhrep(model::JuMP.Model)\n\nBuilds an H-representation from the feasibility set of the JuMP model model. Note that if non-linear constraint are present in the model, they are ignored.\n\n\n\nhrep(A::AbstractMatrix, b::AbstractVector, linset::IntSet=IntSet())\n\nCreates an H-representation for the polyhedron defined by the inequalities langle A_i x rangle = b_i if i in linset and langle A_i x rangle le b_i otherwise where A_i is the ith row of A, i.e. A[i,:] and b_i is b[i].\n\n\n\n"
},

{
    "location": "representation.html#H-representation-1",
    "page": "Representation",
    "title": "H-representation",
    "category": "section",
    "text": "The fundamental element of an H-representation is the halfspaceHalfSpaceAn H-representation can be created as the intersection of several halfspaces. For instance, the polytopebeginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*can be created as follows:HalfSpace([1, 1], 1) ∩ HalfSpace([1, -1], 0) ∩ HalfSpace([-1, 0], 0)Even if HalfSpaces are enough to describe any polyhedron, it is sometimes important to represent the fact that the polyhedron is contained in an affine subspace. For instance, the simplexbeginalign*\n  x_1 + x_2 = 1 \n  x_1 geq 0 \n  x_2 geq 0\nendalign*is 1-dimensional even if it is defined in a 2-dimensional space.The fundamental element of an affine subspace is the hyperplaneHyperPlaneAn affine subspace can be created as the intersection of several hyperplanes. For instanceHyperPlane([1, 1], 1) ∩ HyperPlane([1, 0], 0)represents the 0-dimensional affine subspace only containing the point (0 1).To represent a polyhedron that is not full-dimensional, hyperplanes and halfspaces can be mixed in any order. For instance, the simplex defined above can be obtained as follows:HalfSpace([-1, 0], 0) ∩ HyperPlane([1, 1], 1) ∩ HalfSpace([0, -1], 0)In addition to being created incrementally with intersections, an H-representation can also be created using the hrep functionhrep"
},

{
    "location": "representation.html#Polyhedra.halfspaces",
    "page": "Representation",
    "title": "Polyhedra.halfspaces",
    "category": "function",
    "text": "halfspaces(hrep::HRep)\n\nReturns an iterator over the halfspaces of the H-representation hrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.nhalfspaces",
    "category": "function",
    "text": "nhalfspaces(hrep::HRep)\n\nReturns the number of halfspaces of the H-representation hrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hashalfspaces",
    "page": "Representation",
    "title": "Polyhedra.hashalfspaces",
    "category": "function",
    "text": "hashalfspaces(hrep::HRep)\n\nReturns whether the H-representation hrep has any halfspace.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hyperplanes",
    "page": "Representation",
    "title": "Polyhedra.hyperplanes",
    "category": "function",
    "text": "hyperplanes(hrep::HRep)\n\nReturns an iterator over the hyperplanes of the H-representation hrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nhyperplanes",
    "page": "Representation",
    "title": "Polyhedra.nhyperplanes",
    "category": "function",
    "text": "nhyperplanes(hrep::HRep)\n\nReturns the number of hyperplanes of the H-representation hrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hashyperplanes",
    "page": "Representation",
    "title": "Polyhedra.hashyperplanes",
    "category": "function",
    "text": "hashyperplanes(hrep::HRep)\n\nReturns whether the H-representation hrep has any hyperplane.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.allhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.allhalfspaces",
    "category": "function",
    "text": "allhalfspaces(hrep::HRep)\n\nReturns an iterator over the halfspaces and hyperplanes in the H-representation hrep splitting hyperplanes in two halfspaces.\n\nExamples\n\nhrep = HyperPlane([1, 0], 1) ∩ HalfSpace([0, 1], 1)\ncollect(allhalfspaces(hrep)) # Returns [HalfSpace([1, 0]), HalfSpace([-1, 0]), HalfSpace([0, 1])]\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nallhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.nallhalfspaces",
    "category": "function",
    "text": "nallhalfspaces(hrep::HRep)\n\nReturns the number of halfspaces plus twice the number of hyperplanes in the H-representation hrep, i.e. length(allhalfspaces(hrep))\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasallhalfspaces",
    "page": "Representation",
    "title": "Polyhedra.hasallhalfspaces",
    "category": "function",
    "text": "hasallhalfspaces(hrep::HRep)\n\nReturns whether the H-representation hrep contains any halfspace or hyperplane.\n\n\n\n"
},

{
    "location": "representation.html#Interface-1",
    "page": "Representation",
    "title": "Interface",
    "category": "section",
    "text": "An H-representation is represented as an intersection halfspaces and hyperplanes. The halfspaces can be obtained with halfspaces and the hyperplanes with hyperplanes. As an hyperplane langle a x rangle = beta is the intersection of the two halfspaces langle a x rangle le beta and langle a x rangle ge beta, even if the H-representation contains hyperplanes, a list of halfspaces whose intersection is the polyhedron can be obtained with allhalfspaces, which has nhalfspaces(p) + 2nhyperplanes(p) elements for an H-representation p since each hyperplane is split in two halfspaces.halfspaces\nnhalfspaces\nhashalfspaces\nhyperplanes\nnhyperplanes\nhashyperplanes\nallhalfspaces\nnallhalfspaces\nhasallhalfspaces"
},

{
    "location": "representation.html#Polyhedra.AbstractPoint",
    "page": "Representation",
    "title": "Polyhedra.AbstractPoint",
    "category": "constant",
    "text": "const AbstractPoint{N, T} = Union{Point{N, T}, AbstractVector{T}}\n\nA point in dimension N and of coefficient type T.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.Ray",
    "page": "Representation",
    "title": "Polyhedra.Ray",
    "category": "type",
    "text": "struct Ray{N, T, AT <: MyVec{N, T}}\n    a::AT\nend\n\nThe conic hull of a, i.e. the set of points λa where λ is any nonnegative real number.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.Line",
    "page": "Representation",
    "title": "Polyhedra.Line",
    "category": "type",
    "text": "struct Line{N, T, AT <: MyVec{N, T}}\n    a::AT\nend\n\nThe conic hull of a and -a, i.e. the set of points λa where λ is any real number.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.vrep",
    "page": "Representation",
    "title": "Polyhedra.vrep",
    "category": "function",
    "text": "vrep(p::Polyhedron)\n\nReturns a V-representation for the polyhedron p.\n\n\n\nvrep(lines::LineIt)\n\nCreates an affine space from the list of lines lines.\n\nExamples\n\nvrep([Line([1, 0, 0]), Line([0, 1, 0])])\n\ncreates the 2-dimensional affine subspace containing all the points (x_1 x_2 0), i.e. the x_1x_2-plane.\n\n\n\nvrep(points::PointIt)\n\nCreates a V-representation for the polytope equal to the convex hull of the points points.\n\nExamples\n\nThe convex hull of (0 0), (0 1) and (12 12) can be created as follows using exact arithmetic\n\nvrep([[0, 0], [0, 1], [1//2, 1//2]])\n\nor as follows using floating point arithmetic\n\nvrep([[0, 0], [0, 1], [1/2, 1/2]])\n\n\n\nvrep(lines::LineIt, rays::RayIt)\n\nCreates a V-representation for the polyhedral cone equal to the conic hull of the lines lines and rays rays.\n\nExamples\n\nvrep([Line([0, 1])], [Ray([1, 0])])\n\ncreates a V-representation for the halfspace x_1 ge 0.\n\n\n\nvrep(rays::RayIt)\n\nCreates a V-representation for the polyhedral cone equal to the conic hull of the rays rays.\n\nExamples\n\nvrep([Ray([1, 0]), Ray([0, 1])])\n\ncreates a V-representation for positive orthant.\n\n\n\nvrep(points::PointIt, lines::LineIt, rays::RayIt)\n\nCreates a V-representation for the polyhedron equal to the minkowski sum of the convex hull of points with the conic hull of lines and rays.\n\n\n\nvrep(V::AbstractMatrix, R::AbstractMatrix, Rlinset::IntSet=IntSet())\n\nCreates a V-representation for the polyhedron defined by the points V_i, lines R_i if i in Rlinset and rays R_i otherwise where V_i (resp. R_i) is the ith row of V (resp. R), i.e. V[i,:] (resp. R[i,:]).\n\n\n\n"
},

{
    "location": "representation.html#V-representation-1",
    "page": "Representation",
    "title": "V-representation",
    "category": "section",
    "text": "The fundamental elements of an V-representation are the points and raysAbstractPoint\nRayA V-representation can be created as the minkowski sum between a convex hull of points and a conic hull of rays. For instance, the positive orthant without the simplex defined in the H-representation section can be created as follows:convexhull([1, 0], [0, 1]) + conichull([1, 0], [0, 1])The V-representation represents the polyhedron as a minkowski sum of a polytope and a polyhedral cone. The polytope is represented using a P-representation : a convex hull of points. The polyhedral cone is represented using an R-representation : a conic hull of rays.Even if rays are enough to describe any polyhedral cone, it is sometimes important to represent the fact that the polyhedron contains an affine subspace. For instance, the polyhedron created withconvexhull([1, 0], [0, 1]) + conichull([1, 1], [-1, -1])contains the line [1, 1].The fundamental element of an affine subspace is the lineLineAn affine subspace can be created as the conic hull/minkownski sum of several lines. For instanceconichull(Line([1, 0]), Line([0, 1]))represents the full space.In addition to being created incrementally with convex hull and minkowsky addition, a V-representation can also be created using the vrep functionvrep"
},

{
    "location": "representation.html#Polyhedra.points",
    "page": "Representation",
    "title": "Polyhedra.points",
    "category": "function",
    "text": "points(vrep::VRep)\n\nReturns an iterator over the points of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.npoints",
    "page": "Representation",
    "title": "Polyhedra.npoints",
    "category": "function",
    "text": "npoints(vrep::VRep)\n\nReturns the number of points of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.haspoints",
    "page": "Representation",
    "title": "Polyhedra.haspoints",
    "category": "function",
    "text": "haspoints(vrep::VRep)\n\nReturns whether the V-representation vrep has any point.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.rays",
    "page": "Representation",
    "title": "Polyhedra.rays",
    "category": "function",
    "text": "rays(vrep::VRep)\n\nReturns an iterator over the rays of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nrays",
    "page": "Representation",
    "title": "Polyhedra.nrays",
    "category": "function",
    "text": "nrays(vrep::VRep)\n\nReturns the number of rays of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasrays",
    "page": "Representation",
    "title": "Polyhedra.hasrays",
    "category": "function",
    "text": "hasrays(vrep::VRep)\n\nReturns whether the V-representation vrep has any ray.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.lines",
    "page": "Representation",
    "title": "Polyhedra.lines",
    "category": "function",
    "text": "lines(vrep::VRep)\n\nReturns an iterator over the lines of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nlines",
    "page": "Representation",
    "title": "Polyhedra.nlines",
    "category": "function",
    "text": "nlines(vrep::VRep)\n\nReturns the number of lines of the V-representation vrep.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.haslines",
    "page": "Representation",
    "title": "Polyhedra.haslines",
    "category": "function",
    "text": "haslines(vrep::VRep)\n\nReturns whether the V-representation vrep has any line.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.allrays",
    "page": "Representation",
    "title": "Polyhedra.allrays",
    "category": "function",
    "text": "allrays(vrep::VRep)\n\nReturns an iterator over the rays and lines in the V-representation vrep splitting lines in two rays.\n\nExamples\n\nvrep = Line([1, 0]) + Ray([0, 1])\ncollect(allrays(vrep)) # Returns [Ray([1, 0]), Ray([-1, 0]), Ray([0, 1])]\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nallrays",
    "page": "Representation",
    "title": "Polyhedra.nallrays",
    "category": "function",
    "text": "nallrays(vrep::VRep)\n\nReturns the number of rays plus twice the number of lines in the V-representation vrep, i.e. length(allrays(vrep))\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasallrays",
    "page": "Representation",
    "title": "Polyhedra.hasallrays",
    "category": "function",
    "text": "hasallrays(vrep::VRep)\n\nReturns whether the V-representation vrep contains any ray or line.\n\n\n\n"
},

{
    "location": "representation.html#Interface-2",
    "page": "Representation",
    "title": "Interface",
    "category": "section",
    "text": "A P-representation is represented as a convex hull points. The points can be obtained with points.points\nnpoints\nhaspointsAn R-representation is represented as a conic hull of lines and rays. The rays can be obtained with rays and the lines with lines. As a line r is the conic hull of of the two rays r and -r, even if the R-representation contains lines, a list of rays whose conic hull is the polyhedral cone can be obtained with allrays, which has nrays(R) + 2nlines(R) elements for an R-representation R since each line is split in two rays.rays\nnrays\nhasrays\nlines\nnlines\nhaslines\nallrays\nnallrays\nhasallrays"
},

{
    "location": "polyhedron.html#",
    "page": "Polyhedron",
    "title": "Polyhedron",
    "category": "page",
    "text": ""
},

{
    "location": "polyhedron.html#Polyhedron-1",
    "page": "Polyhedron",
    "title": "Polyhedron",
    "category": "section",
    "text": "As seen in the previous section, a polyhedron can be described in 2 ways: either using the H-representation (list of inequalities) or the V-representation (list of points and rays). A typical problem is: Given the H-(or V-)representation of one or several polyhedra, what is the H-(or V-)representation of some polyhedra obtained after some operations on these initial polyhedra. This description is similar to the description usually given to algorithms except that in that case we talk about numbers given in their binary representation and not polyhedra given in their H-(or V-)representation. This motivates the creation of a type representing polyhedra. Just like the abstract type AbstractArray{N,T} represents an N-dimensional array with elements of type T, the abstract type Polyhedron{N,T} represents an N-dimensional polyhedron with elements of coefficient type T.There is typically one concrete subtype of Polyhedron by library. For instance, the CDD library defines CDDPolyhedron and the LRS library defines LRSPolyhedron. It must be said that the type T is not necessarily how the elements are stored internally by the library but the polyhedron will behave just like it is stored that way. For instance, when retreiving an H-(or V-)representation, the representation will be of type T. Therefore using Int for T may result in InexactError. For this reason, by default, the type T chosen is not a subtype of Integer.Consider the representations hrep, vrep and vrepf created in the preceding section. One can use the CDD library, to create an instance of a concrete subtype of Polyhedronjulia> using CDDLib\njulia> polyf = polyhedron(hrep, CDDLibrary())\njulia> typeof(polyhf)\nCDDLib.CDDPolyhedron{2,Float64}We see that the library has choosen to deal with floating point arithmetic. This decision does not depend on the type of hrep but only on the instance of CDDLibrary given. CDDLibrary creates CDDPolyhedron of type either Float64 or Rational{BigInt}. One can choose the first one using CDDLibrary(:float) and the second one using CDDLibrary(:exact), by default it is :float.julia> poly = polyhedron(hrep, CDDLibrary(:exact))\njulia> typeof(poly)\nCDDLib.CDDPolyhedron{2,Rational{BigInt}}The first polyhedron polyf can also be created from its V-representation using either of the 4 following linesjulia> polyf = polyhedron(vrepf, CDDLibrary(:float))\njulia> polyf = polyhedron(vrepf, CDDLibrary())\njulia> polyf = polyhedron(vrep,  CDDLibrary(:float))\njulia> polyf = polyhedron(vrep,  CDDLibrary())and poly using either of those linesjulia> poly = polyhedron(vrepf, CDDLibrary(:exact))\njulia> poly = polyhedron(vrep , CDDLibrary(:exact))of course, creating a representation in floating points with exact arithmetic works here because we have 0.5 which is 0.1 in binary but in general, is not a good idea.julia> Rational{BigInt}(1/2)\n1//2\njulia> Rational{BigInt}(1/3)\n6004799503160661//18014398509481984\njulia> Rational{BigInt}(1/5)\n3602879701896397//18014398509481984"
},

{
    "location": "polyhedron.html#Retrieving-a-representation-1",
    "page": "Polyhedron",
    "title": "Retrieving a representation",
    "category": "section",
    "text": "One can retrieve an H-representation (resp. V-representation) from a polyhedron using hrep (resp. vrep). The concrete subtype of HRepresentation (resp. VRepresentation) returned is not necessarily the same that the one used to create the polyhedron. As a rule of thumb, it is the representation the closest to the internal representation used by the library.julia> hrep = hrep(poly)\njulia> typeof(hrep)\nPolyhedra.LiftedHRepresentation{2,Rational{BigInt}}\njulia> hrep = SimpleHRepresentation(hrep)\njulia> typeof(hrep)\nPolyhedra.SimpleHRepresentation{2,Rational{BigInt}}\njulia> hrep.A\n3x2 Array{Rational{BigInt},2}:\n  1//1   1//1\n  1//1  -1//1\n -1//1   0//1\njulia> hrep.b\n3-element Array{Rational{BigInt},1}:\n 1//1\n 0//1\n 0//1\njulia> vrep = vrep(poly)\njulia> typeof(vrep)\nPolyhedra.LiftedVRepresentation{2,Rational{BigInt}}\njulia> vrep = SimpleVRepresentation(vrep)\njulia> typeof(vrep)\nPolyhedra.SimpleVRepresentation{2,Rational{BigInt}}\njulia> vrep.V\n3x2 Array{Rational{BigInt},2}:\n 1//2  1//2\n 0//1  1//1\n 0//1  0//1\n\njulia> vrep.R\n0x2 Array{Rational{BigInt},2}"
},

{
    "location": "polyhedron.html#Creating-a-polyhedron-from-the-feasible-set-of-a-JuMP-model-1",
    "page": "Polyhedron",
    "title": "Creating a polyhedron from the feasible set of a JuMP model",
    "category": "section",
    "text": "A typical application of polyhedral computation is the computation of the set of extreme points and rays of the feasible set of an optimization problem. This comes from the fact that given a minimization of a concave function (or maximization of a convex function) on a convex feasible set (e.g. Linear Programming), we are either in the following three situations:The feasible set is empty, i.e. the problem is infeasible.\nAn extreme ray is optimal, i.e. the problem is unbounded (or it may also be bounded if the objective is constant along the ray).\nAn extreme point is optimal.A JuMP model is treated by polyhedron just like any H-representation. For example, the hypercube of dimension n can be created as followsm = Model()\n@variable(m, 0 ≤ x[1:n] ≤ 1)\n\npoly = polyhedron(m, CDDLibrary(:exact))In fact, the MathProgBase representation of the feasible set of a linear program:beginalign*\n  lb leq Ax leq ub\n  l leq x leq u\nendalign*has LPHRepresentation as a corresponding H-representation. A JuMP Model can be converted to this representation using LPHRepresentation(m)."
},

{
    "location": "polyhedron.html#Polyhedra.hrepiscomputed",
    "page": "Polyhedron",
    "title": "Polyhedra.hrepiscomputed",
    "category": "function",
    "text": "hrepiscomputed(p::Polyhedron)\n\nReturns whether the H-representation of this polyhedron has been computed.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.vrepiscomputed",
    "page": "Polyhedron",
    "title": "Polyhedra.vrepiscomputed",
    "category": "function",
    "text": "vrepiscomputed(p::Polyhedron)\n\nReturns whether the V-representation of this polyhedron has been computed.\n\n\n\n"
},

{
    "location": "polyhedron.html#Checking-if-a-representation-has-been-computed-1",
    "page": "Polyhedron",
    "title": "Checking if a representation has been computed",
    "category": "section",
    "text": "hrepiscomputed\nvrepiscomputed"
},

{
    "location": "polyhedron.html#Polyhedra.incidenthalfspaces",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthalfspaces",
    "category": "function",
    "text": "incidenthalfspaces(p::Polyhedron, idx)\n\nReturns the list of halfspaces incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidenthalfspaceindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthalfspaceindices",
    "category": "function",
    "text": "incidenthalfspaceindices(p::Polyhedron, idx)\n\nReturns the list of the indices of halfspaces incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentpoints",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentpoints",
    "category": "function",
    "text": "incidentpoints(p::Polyhedron, idx)\n\nReturns the list of points incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentpointindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentpointindices",
    "category": "function",
    "text": "incidentpointindices(p::Polyhedron, idx)\n\nReturns the list of the indices of points incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentrays",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentrays",
    "category": "function",
    "text": "incidentrays(p::Polyhedron, idx)\n\nReturns the list of rays incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentrayindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentrayindices",
    "category": "function",
    "text": "incidentrayindices(p::Polyhedron, idx)\n\nReturns the list of the indices of rays incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidenthyperplanes",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthyperplanes",
    "category": "function",
    "text": "incidenthyperplanes(p::Polyhedron, idx)\n\nReturns the list of hyperplanes incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidenthyperplaneindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidenthyperplaneindices",
    "category": "function",
    "text": "incidenthyperplaneindices(p::Polyhedron, idx)\n\nReturns the list of the indices of hyperplanes incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentlines",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentlines",
    "category": "function",
    "text": "incidentlines(p::Polyhedron, idx)\n\nReturns the list of lines incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Polyhedra.incidentlineindices",
    "page": "Polyhedron",
    "title": "Polyhedra.incidentlineindices",
    "category": "function",
    "text": "incidentlineindices(p::Polyhedron, idx)\n\nReturns the list of the indices of lines incident to idx for the polyhedron p.\n\n\n\n"
},

{
    "location": "polyhedron.html#Incidence-1",
    "page": "Polyhedron",
    "title": "Incidence",
    "category": "section",
    "text": "A point p (ray r) is incident to an halfspace langle a x rangle le beta if langle a p rangle = beta (resp. langle a r rangle = beta).incidenthalfspaces\nincidenthalfspaceindices\nincidentpoints\nincidentpointindices\nincidentrays\nincidentrayindicesIn a polyhedron, all points and rays are incident to all hyperplanes and all halfspaces are incident to all lines. The following methods are therefore redundant, e.g. incidenthyperplanes(p, idx) is equivalent to hyperplanes(p) and incidenthyperplaneindices(p, idx) is equivalent to eachindex(hyperplanes(p)). The methods are hence only defined for consistency.incidenthyperplanes\nincidenthyperplaneindices\nincidentlines\nincidentlineindices"
},

{
    "location": "redundancy.html#",
    "page": "Containment/Redundancy",
    "title": "Containment/Redundancy",
    "category": "page",
    "text": ""
},

{
    "location": "redundancy.html#Containment/Redundancy-1",
    "page": "Containment/Redundancy",
    "title": "Containment/Redundancy",
    "category": "section",
    "text": ""
},

{
    "location": "redundancy.html#Base.in",
    "page": "Containment/Redundancy",
    "title": "Base.in",
    "category": "function",
    "text": "in(p::VRepElement, h::HRepElement)\n\nReturns whether p is in h. If h is an hyperplane, it returns whether langle a x rangle approx beta. If h is an halfspace, it returns whether langle a x rangle le beta.\n\nin(p::VRepElement, h::HRep)\n\nReturns whether p is in h, e.g. in all the hyperplanes and halfspaces supporting h.\n\n\n\n"
},

{
    "location": "redundancy.html#Base.issubset",
    "page": "Containment/Redundancy",
    "title": "Base.issubset",
    "category": "function",
    "text": "issubset(p::Rep, h::HRepElement)\n\nReturns whether p is a subset of h, i.e. whether h supports the polyhedron p.\n\n\n\n"
},

{
    "location": "redundancy.html#Polyhedra.ininterior",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.ininterior",
    "category": "function",
    "text": "ininterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the interior of h. If h is an hyperplane, it always returns false. If h is an halfspace langle a x rangle leq beta, it returns whether p is in the open halfspace langle a x rangle  beta\n\nininterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the interior of h, e.g. in the interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n"
},

{
    "location": "redundancy.html#Polyhedra.inrelativeinterior",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.inrelativeinterior",
    "category": "function",
    "text": "inrelativeinterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the relative interior of h. If h is an hyperplane, it is equivalent to p in h since the relative interior of an hyperplane is itself. If h is an halfspace, it is equivalent to ininterior(p, h).\n\ninrelativeinterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the relative interior of h, e.g. in the relative interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n"
},

{
    "location": "redundancy.html#Containment-1",
    "page": "Containment/Redundancy",
    "title": "Containment",
    "category": "section",
    "text": "in\nissubset\nininterior\ninrelativeinterior"
},

{
    "location": "redundancy.html#Polyhedra.dim",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.dim",
    "category": "function",
    "text": "dim(p::Polyhedron)\n\nReturns the dimension of the affine hull of the polyhedron. That is the number of non-redundant hyperplanes that define it.\n\n\n\n"
},

{
    "location": "redundancy.html#Polyhedra.isredundant",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.isredundant",
    "category": "function",
    "text": "isredundant(p::Rep, idx::Index; strongly=false)\n\nReturn a Bool indicating whether the element with index idx can be removed without changing the polyhedron represented by p. If strongly is true,\n\nif idx is an H-representation element h, it returns true only if no V-representation element of p is in the hyperplane of h.\nif idx is a V-representation element v, it returns true only if v is in the relative interior of p.\n\n\n\n"
},

{
    "location": "redundancy.html#Polyhedra.removehredundancy!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removehredundancy!",
    "category": "function",
    "text": "removehredundancy!(p::HRep)\n\nRemoves the elements of the H-representation of p that can be removed without changing the polyhedron represented by p. That is, it only keeps the halfspaces corresponding to facets of the polyhedron.\n\n\n\n"
},

{
    "location": "redundancy.html#Polyhedra.removevredundancy!",
    "page": "Containment/Redundancy",
    "title": "Polyhedra.removevredundancy!",
    "category": "function",
    "text": "removevredundancy!(p::VRep)\n\nRemoves the elements of the V-representation of p that can be removed without changing the polyhedron represented by p. That is, it only keeps the extreme points and rays. This operation is often called \"convex hull\" as the remaining points are the extreme points of the convex hull of the initial set of points.\n\n\n\n"
},

{
    "location": "redundancy.html#Redundancy-1",
    "page": "Containment/Redundancy",
    "title": "Redundancy",
    "category": "section",
    "text": "dim\nisredundant\nremovehredundancy!\nremovevredundancy!"
},

{
    "location": "projection.html#",
    "page": "Projection/Elimination",
    "title": "Projection/Elimination",
    "category": "page",
    "text": ""
},

{
    "location": "projection.html#Polyhedra.FourierMotzkin",
    "page": "Projection/Elimination",
    "title": "Polyhedra.FourierMotzkin",
    "category": "type",
    "text": "FourierMotzkin\n\nComputation of the projection by computing the H-representation and applying the Fourier-Motzkin elimination algorithm to it.\n\n\n\n"
},

{
    "location": "projection.html#Polyhedra.BlockElimination",
    "page": "Projection/Elimination",
    "title": "Polyhedra.BlockElimination",
    "category": "type",
    "text": "BlockElimination\n\nComputation of the projection by computing the H-representation and applying the block elimination algorithm to it.\n\n\n\n"
},

{
    "location": "projection.html#Polyhedra.ProjectGenerators",
    "page": "Projection/Elimination",
    "title": "Polyhedra.ProjectGenerators",
    "category": "type",
    "text": "ProjectGenerators\n\nComputation of the projection by computing the V-representation and projecting them.\n\n\n\n"
},

{
    "location": "projection.html#Polyhedra.eliminate",
    "page": "Projection/Elimination",
    "title": "Polyhedra.eliminate",
    "category": "function",
    "text": "eliminate(p::Polyhedron, delset, algo::EliminationAlgorithm)\n\nEliminate the dimensions in delset by projecting the polyhedron onto the remaining dimension.\n\n\n\n"
},

{
    "location": "projection.html#Polyhedra.project",
    "page": "Projection/Elimination",
    "title": "Polyhedra.project",
    "category": "function",
    "text": "project(p::Polyhedron, pset, algo)\n\nEquivalent to `eliminate(p, setdiff(1:N, pset), algo).\n\n\n\n"
},

{
    "location": "projection.html#Polyhedra.fixandeliminate",
    "page": "Projection/Elimination",
    "title": "Polyhedra.fixandeliminate",
    "category": "function",
    "text": "fixandeliminate(p::HRep{N, T}, I, v)\n\nFix the variables with indices in I to the corresponding value in v. This is equivalent to doing the following:\n\nfunction ei(i)\n    a = zeros(T, N)\n    a[i] = one(T)\n    a\nend\neliminate(p ∩ HyperPlane(ei(I[1], v[1]) ∩ ... ∩ HyperPlane(ei(I[1], v[1]))\n\nbut it is much more efficient. The code above does a polyhedral projection while this function simply replace each halfspace ⟨a, x⟩ ≤ β (resp. each hyperplane ⟨a, x⟩ = β) by the halfspace ⟨a_J, x⟩ ≤ β - ⟨a_I, v⟩ (resp. the hyperplane ⟨a_J, x⟩ = β - ⟨a_I, v⟩) where J = setdiff(1:N, I).\n\n\n\n"
},

{
    "location": "projection.html#Projection/Elimination-1",
    "page": "Projection/Elimination",
    "title": "Projection/Elimination",
    "category": "section",
    "text": "Consider the polyhedron created in the beginning of this section. As a reminder, it represents the following H-representation:beginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*One can verify that for any 0 leq x_2 leq 1, there exists a value x_1 such that (x_1 x_2) is in this polyhedron. This means that the H-representation obtained by eliminating x_1 is:beginalign*\n  x_1  leq 1 \n  x_1  geq 0\nendalign*where x_1 in the H-representation above represents x_2 in the previous one. This can be obtained as followsjulia> poly_x2 = eliminate(poly, [1])\njulia> hrep(poly_x2)\nH-representation\nbegin\n 2 2 rational\n 1//1 -1//1\n 0//1 1//1\nendThere is two methods of computing the elimination implemented in CDDLib: Fourier-Motzkin elimination and block elimination. As written by K. Fukuda in CDD\'s documentation, \"[Block elimination] might be a faster way to eliminate variables than the repeated [Fourier-Motzkin elimination] when the number of variables to eliminate is large\". You can specify the method to use as a third argument, e.g. eliminate(poly, [1], FourierMotzkin()), eliminate(poly, [1], BlockElimination()). A third method can be chosen: ProjectGenerators. It computes the V-representation and then project each of its elements. This is the method of choice when the V-representation is already computed.FourierMotzkin\nBlockElimination\nProjectGeneratorsIf nothing is specified as in the block of code above, the behavior depends on the polyhedral library. If neither Fourier-Motzkin nor block elimination is implemented or if the V-representation is already computed then :ProjectGenerators is chosen. Otherwise, Polyhedra lets the library decide. In CDDLib, :FourierMotzkin is chosen when only the last dimension needs to be eliminated and :BlockElimination is chosen otherwise. Note that CDDLib only supports projecting the last trailing dimensions.eliminate\nproject\nfixandeliminate"
},

{
    "location": "utilities.html#",
    "page": "Utilities",
    "title": "Utilities",
    "category": "page",
    "text": ""
},

{
    "location": "utilities.html#Utilities-1",
    "page": "Utilities",
    "title": "Utilities",
    "category": "section",
    "text": ""
},

{
    "location": "utilities.html#Base.:*",
    "page": "Utilities",
    "title": "Base.:*",
    "category": "function",
    "text": "*(p1::Rep, p2::Rep)\n\nCartesian product between the polyhedra p1 and p2.\n\n\n\n*(P::AbstractMatrix, p::VRep)\n\nTransform the polyhedron represented by p into P p by transforming each element of the V-representation (points, symmetric points, rays and lines) x into P x.\n\n\n\n"
},

{
    "location": "utilities.html#Base.:\\",
    "page": "Utilities",
    "title": "Base.:\\",
    "category": "function",
    "text": "\\(P::AbstractMatrix, p::HRep)\n\nTransform the polyhedron represented by p into P^-1 p by transforming each halfspace langle a x rangle le beta into langle P^top a x rangle le beta and each hyperplane langle a x rangle = beta into langle P^top a x rangle = beta.\n\n\n\n"
},

{
    "location": "utilities.html#Base.:/",
    "page": "Utilities",
    "title": "Base.:/",
    "category": "function",
    "text": "/(p::HRep, P::AbstractMatrix)\n\nTransform the polyhedron represented by p into P^-T p by transforming each halfspace langle a x rangle le beta into langle P a x rangle le beta and each hyperplane langle a x rangle = beta into langle P a x rangle = beta.\n\n\n\n"
},

{
    "location": "utilities.html#Base.intersect",
    "page": "Utilities",
    "title": "Base.intersect",
    "category": "function",
    "text": "intersect(P1::HRep, P2::HRep)\n\nTakes the intersection of P1 and P2  x  x in P_1 x in P_2 . It is very efficient between two H-representations or between two polyhedron for which the H-representation has already been computed. However, if P1 (resp. P2) is a polyhedron for which the H-representation has not been computed yet, it will trigger a representation conversion which is costly. See the Polyhedral Computation FAQ for a discussion on this operation.\n\nThe type of the result will be chosen closer to the type of P1. For instance, if P1 is a polyhedron (resp. H-representation) and P2 is a H-representation (resp. polyhedron), intersect(P1, P2) will be a polyhedron (resp. H-representation). If P1 and P2 are both polyhedra (resp. H-representation), the resulting polyhedron type (resp. H-representation type) will be computed according to the type of P1. The coefficient type however, will be promoted as required taking both the coefficient type of P1 and P2 into account.\n\n\n\nintersect(v::VRep{N, T}, h::HRepElement)\n\nCompute the intersection of v with an halfspace or hyperplane h. The method used by default is to keep the V-representation element of v that are in h and add new ones generated as the intersection between the hyperplane defining h and the segment between two adjacent V-representation elements of v that are in either sides of the hyperplane. See Lemma 3 of [1] for more detail on the method.\n\n[1] Fukuda, K. and Prodon, A. Double description method revisited Combinatorics and computer science, Springer, 1996, 91-111\n\n\n\n"
},

{
    "location": "utilities.html#Base.intersect!",
    "page": "Utilities",
    "title": "Base.intersect!",
    "category": "function",
    "text": "intersect!(p1::VRep, p2::VRep)\n\nSame as intersect except that p1 is modified to be equal to the intersection.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.convexhull",
    "page": "Utilities",
    "title": "Polyhedra.convexhull",
    "category": "function",
    "text": "convexhull(P1::VRep, P2::VRep)\n\nTakes the convex hull of P1 and P2  lambda x + (1-lambda) y  x in P_1 y in P_2 . It is very efficient between two V-representations or between two polyhedron for which the V-representation has already been computed. However, if P1 (resp. P2) is a polyhedron for which the V-representation has not been computed yet, it will trigger a representation conversion which is costly.\n\nThe type of the result will be chosen closer to the type of P1. For instance, if P1 is a polyhedron (resp. V-representation) and P2 is a V-representation (resp. polyhedron), convexhull(P1, P2) will be a polyhedron (resp. V-representation). If P1 and P2 are both polyhedra (resp. V-representation), the resulting polyhedron type (resp. V-representation type) will be computed according to the type of P1. The coefficient type however, will be promoted as required taking both the coefficient type of P1 and P2 into account.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.convexhull!",
    "page": "Utilities",
    "title": "Polyhedra.convexhull!",
    "category": "function",
    "text": "convexhull!(p1::VRep, p2::VRep)\n\nSame as convexhull except that p1 is modified to be equal to the convex hull.\n\n\n\n"
},

{
    "location": "utilities.html#Operations-1",
    "page": "Utilities",
    "title": "Operations",
    "category": "section",
    "text": "*\n\\\n/\nintersect\nintersect!\nconvexhull\nconvexhull!"
},

{
    "location": "utilities.html#Polyhedra.volume",
    "page": "Utilities",
    "title": "Polyhedra.volume",
    "category": "function",
    "text": "volume(p::Polyhedron{N, T}) where {N, T}\n\nReturns the N-dimensional hyper-volume of the polyhedron p. Returns Inf or -one(T) if it is infinite depending on whether the type T has an infinite value.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.surface",
    "page": "Utilities",
    "title": "Polyhedra.surface",
    "category": "function",
    "text": "surface(p::Polyhedron{N, T}) where {N, T}\n\nReturns the N-1-dimensional hyper-volume of the surface of the polyhedron p. Returns Inf or -one(T) if it is infinite depending on whether the type T has an infinite value.\n\n\n\n"
},

{
    "location": "utilities.html#Volume-1",
    "page": "Utilities",
    "title": "Volume",
    "category": "section",
    "text": "volume\nsurface"
},

{
    "location": "utilities.html#Polyhedra.chebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.chebyshevcenter",
    "category": "function",
    "text": "chebyshevcenter(p::Rep[, solver])\n\nIf p is a H-representation or is a polyhedron for which the H-representation has already been computed, calls hchebyshevcenter, otherwise, call vchebyshevcenter.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.hchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.hchebyshevcenter",
    "category": "function",
    "text": "hchebyshevcenter(p::HRep[, solver])\n\nReturn a tuple with the center and radius of the largest euclidean ball contained in the polyhedron p. Throws an error if the polyhedron is empty or if the radius is infinite.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.vchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.vchebyshevcenter",
    "category": "function",
    "text": "vchebyshevcenter(p::VRep[, solver])\n\nReturn a tuple with the center and radius of the smallest euclidean ball containing the polyhedron p. Throws an error if the polyhedron is empty or if the radius is infinite (i.e. p is not a polytope, it contains rays).\n\n\n\n"
},

{
    "location": "utilities.html#Chebyshev-center-1",
    "page": "Utilities",
    "title": "Chebyshev center",
    "category": "section",
    "text": "chebyshevcenter\nhchebyshevcenter\nvchebyshevcenter"
},

]}
