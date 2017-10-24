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
    "text": "Polyhedra is a package for polyhedra manipulations in Julia. It provides an unified interface for Polyhedra Manipulation Libraries such as CDDLib.jl and LRSLib.jl.Polyhedra can either be represented by a set of linear inequalities or by vertices and rays. In the first case, the points of the polyhedron are the points which satisfies all the inequalities and in the second case they are the points that can be expressed as a convex combination of the vertices plus a conic combination of the rays. The manipulations that Polyhedra can perform includeProjection: Projection of a polyhedron on a lower dimensional space, e.g. Fourier-Motzkin elimination.\nChanging the Representation\nVertex enumeration problem: Computing the extremal vertices and rays from an inequality representation\nConvex hull problem: Computing a set of linear inequalities describing the polyhedron from a vertex/ray representation\nRemoval of redundant inequalities or redundant vertices/rays.\nDecomposition of 3D/2D polyhedra into a points and triangular faces, enabling easy visualization of 3D/2D polyhedra using DrakeVisualizer or GLVisualize.Depending on the library, those manipulation can either be in floating point or exact rational arithmetic.Polyhedra remains under active development, and we welcome your feedback, suggestions, and bug reports."
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
    "text": "Pages = [\"installation.md\", \"representation.md\", \"polyhedron.md\", \"utilities.md\"]\nDepth = 2"
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
    "text": "Open a Julia console (e.g. enter julia at the command line) and writejulia> Pkg.add(\"Polyhedra\")To start using Polyhedra, you can now just writejulia> using PolyhedraPolyhedra includes a default implementation for every operation but external library can also be used. See the next section on installing a library."
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
    "location": "representation.html#Representation-1",
    "page": "Representation",
    "title": "Representation",
    "category": "section",
    "text": "Polyhedra can be described in 2 different ways.H-representation: As the intersection of finitely many halfspaces given by its facets.\nV-representation: As the convex hull of its vertices + the conic hull of its rays where '+' is the Minkowski sum.In Polyhedra.jl, those representations are given the respective abstract types HRepresentation and VRepresentation which are themself subtypes of Representation.For instance, consider the 2-dimensional polyhedron described by the following H-representation:beginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*This set of inequalities can be written in the matrix form Ax leq b whereA = beginpmatrix1  11  -1-1  0endpmatrix b = beginpmatrix100endpmatrixLet's create this H-representation using the concrete subtype SimpleHRepresentation of the abstract type HRepresentation.julia> using Polyhedra\njulia> A = [1 1;1 -1;-1 0]\njulia> b = [1,0,0]\njulia> hrep = SimpleHRepresentation(A, b)\njulia> typeof(hrep)\nPolyhedra.SimpleHRepresentation{2,Int64}This polyhedron has three vertices: (00), (01) and (0505). We can create this V-representation using the concrete subtype SimpleVRepresentation of the abstract type VRepresentation. Because 05 is fractional, have two choices: either use exact rational arithemticjulia> V = [0 0; 0 1; 1//2 1//2]\njulia> vrep = SimpleVRepresentation(V)\njulia> typeof(vrep)\nPolyhedra.SimpleVRepresentation{2,Rational{Int64}}or use floating point arithmeticjulia> Vf = [0 0; 0 1; 1/2 1/2]\njulia> vrepf = SimpleVRepresentation(Vf)\njulia> typeof(vrepf)\nPolyhedra.SimpleVRepresentation{2,Float64}"
},

{
    "location": "representation.html#Polyhedra.fulldim",
    "page": "Representation",
    "title": "Polyhedra.fulldim",
    "category": "Function",
    "text": "fulldim(rep::Rep)\n\nReturns the dimension of the space in which the representation is defined. That is, a straight line in a 3D space has fulldim 3.\n\n\n\n"
},

{
    "location": "representation.html#Representation-interface-1",
    "page": "Representation",
    "title": "Representation interface",
    "category": "section",
    "text": "These functions can be called on both H-representation and V-representationfulldim"
},

{
    "location": "representation.html#Polyhedra.hreps",
    "page": "Representation",
    "title": "Polyhedra.hreps",
    "category": "Function",
    "text": "hreps(hr::HRep)\n\nReturns an iterator over the elements of the H-representation.\n\nNote\n\nThis is type unstable as the iterator returns both halfspaces and hyperplanes. It is therefore more efficient to call eqs and ineqs separately.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nhreps",
    "page": "Representation",
    "title": "Polyhedra.nhreps",
    "category": "Function",
    "text": "nhreps(hr::HRep)\n\nReturns the number of halfspaces and hyperplanes of the H-representation.\n\nNote\n\nNote that it does not do redundancy removal so it is not the minimal number of halfspace and hyperplanes needed to represent the polyhedron, it is simply the number that are currently used.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hashreps",
    "page": "Representation",
    "title": "Polyhedra.hashreps",
    "category": "Function",
    "text": "hashreps(hr::HRep)\n\nReturns whether the H-representation contain any halfspace or hyperplane.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.eqs",
    "page": "Representation",
    "title": "Polyhedra.eqs",
    "category": "Function",
    "text": "eqs(hr::HRep)\n\nReturns an iterator over the hyperplanes of the H-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.neqs",
    "page": "Representation",
    "title": "Polyhedra.neqs",
    "category": "Function",
    "text": "neqs(hr::HRep)\n\nReturns the number of hyperplanes of the H-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.haseqs",
    "page": "Representation",
    "title": "Polyhedra.haseqs",
    "category": "Function",
    "text": "haseqs(hr::HRep)\n\nReturns whether the H-representation contain any hyperplane.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.ineqs",
    "page": "Representation",
    "title": "Polyhedra.ineqs",
    "category": "Function",
    "text": "ineqs(hr::HRep)\n\nReturns an iterator over the halfspaces of the H-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nineqs",
    "page": "Representation",
    "title": "Polyhedra.nineqs",
    "category": "Function",
    "text": "nineqs(hr::HRep)\n\nReturns the number of halfspaces of the H-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasineqs",
    "page": "Representation",
    "title": "Polyhedra.hasineqs",
    "category": "Function",
    "text": "hasineqs(hr::HRep)\n\nReturns whether the H-representation contain any halfspace.\n\n\n\n"
},

{
    "location": "representation.html#H-representation-interface-1",
    "page": "Representation",
    "title": "H-representation interface",
    "category": "section",
    "text": "hreps\nnhreps\nhashreps\neqs\nneqs\nhaseqs\nineqs\nnineqs\nhasineqs"
},

{
    "location": "representation.html#Polyhedra.vreps",
    "page": "Representation",
    "title": "Polyhedra.vreps",
    "category": "Function",
    "text": "vreps(vr::VRep)\n\nReturns an iterator over the elements of the V-representation.\n\nNote\n\nThis is type unstable as the iterator returns both points and rays. It is therefore more efficient to call points and rays separately.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nvreps",
    "page": "Representation",
    "title": "Polyhedra.nvreps",
    "category": "Function",
    "text": "nvreps(vr::VRep)\n\nReturns the number of points and rays of the V-representation.\n\nNote\n\nNote that it does not do redundancy removal so it is not the minimal number of points and rays needed to represent the polyhedron, it is simply the number that are currently used.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasvreps",
    "page": "Representation",
    "title": "Polyhedra.hasvreps",
    "category": "Function",
    "text": "hasvreps(vr::VRep)\n\nReturns whether the V-representation contain any ray or point.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.rays",
    "page": "Representation",
    "title": "Polyhedra.rays",
    "category": "Function",
    "text": "rays(vr::VRep)\n\nReturns an iterator over the rays of the V-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.nrays",
    "page": "Representation",
    "title": "Polyhedra.nrays",
    "category": "Function",
    "text": "nrays(hr::HRep)\n\nReturns the number of rays of the V-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.hasrays",
    "page": "Representation",
    "title": "Polyhedra.hasrays",
    "category": "Function",
    "text": "hasrays(vr::VRep)\n\nReturns whether the V-representation contain any ray.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.points",
    "page": "Representation",
    "title": "Polyhedra.points",
    "category": "Function",
    "text": "points(vr::VRep)\n\nReturns an iterator over the points of the V-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.npoints",
    "page": "Representation",
    "title": "Polyhedra.npoints",
    "category": "Function",
    "text": "npoints(vr::VRep)\n\nReturns the number of points of the V-representation.\n\n\n\n"
},

{
    "location": "representation.html#Polyhedra.haspoints",
    "page": "Representation",
    "title": "Polyhedra.haspoints",
    "category": "Function",
    "text": "haspoints(vr::VRep)\n\nReturns whether the V-representation contain any point.\n\n\n\n"
},

{
    "location": "representation.html#V-representation-interface-1",
    "page": "Representation",
    "title": "V-representation interface",
    "category": "section",
    "text": "vreps\nnvreps\nhasvreps\nrays\nnrays\nhasrays\npoints\nnpoints\nhaspoints"
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
    "text": "As seen in the previous section, a polyhedron can be described in 2 ways: either using the H-representation or the V-representation. A typical problem is: Given the H-(or V-)representation of one or several polyhedra, what is the H-(or V-)representation of some polyhedra obtained after some operations of those initial polyhedra. This description is similar to the description usually given to algorithms except that in that case we talk about numbers given in their binary representation and not polyhedra given in their H-(or V-)representation. This motivates the creation of a type representing polyhedra. Just like the abstract type AbstractArray{N,T} represents an N-dimensional array of elements of type T, the abstract type Polyhedron{N,T} represents an N-dimensional polyhedron of elements of type T.There is typically one concrete subtype of Polyhedron by library. For instance, the CDD library defines CDDPolyhedron and the LRS library defines LRSPolyhedron. It must be said that the type T is not necessarily how the elements are stored internally by the library but the polyhedron will behave just like it is stored that way. For instance, when retreiving an H-(or V-)representation, the representation will be of type T. Therefore Int for T is may result in InexactError. For this reason, by default, the type T chosen is not a subtype of Integer.Consider the representations hrep, vrep and vrepf created in the preceding section. One can use the CDD library, to create an instance of a concrete subtype of Polyhedronjulia> using CDDLib\njulia> polyf = polyhedron(hrep, CDDLibrary())\njulia> typeof(polyhf)\nCDDLib.CDDPolyhedron{2,Float64}We see that the library has choosen to deal with floating point arithmetic. This decision does not depend on the type of hrep but only on the instance of CDDLibrary given. CDDLibrary creates CDDPolyhedron of type either Float64 or Rational{BigInt}. One can choose the first one using CDDLibrary(:float) and the second one using CDDLibrary(:exact), by default it is :float.julia> poly = polyhedron(hrep, CDDLibrary(:exact))\njulia> typeof(poly)\nCDDLib.CDDPolyhedron{2,Rational{BigInt}}The first polyhedron polyf can also be created from its V-representation using either of the 4 following linesjulia> polyf = polyhedron(vrepf, CDDLibrary(:float))\njulia> polyf = polyhedron(vrepf, CDDLibrary())\njulia> polyf = polyhedron(vrep,  CDDLibrary(:float))\njulia> polyf = polyhedron(vrep,  CDDLibrary())and poly using either of those linesjulia> poly = polyhedron(vrepf, CDDLibrary(:exact))\njulia> poly = polyhedron(vrep , CDDLibrary(:exact))of course, creating a representation in floating points with exact arithmetic works here because we have 0.5 which is 0.1 in binary but in general, is not a good idea.julia> Rational{BigInt}(1/2)\n1//2\njulia> Rational{BigInt}(1/3)\n6004799503160661//18014398509481984\njulia> Rational{BigInt}(1/5)\n3602879701896397//18014398509481984"
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
    "location": "polyhedron.html#Projecting-a-polyhedron-1",
    "page": "Polyhedron",
    "title": "Projecting a polyhedron",
    "category": "section",
    "text": "Consider the polyhedron created in the beginning of this section. As a reminder, it represents the following H-representation:beginalign*\n  x_1 + x_2 leq 1 \n  x_1 - x_2 leq 0 \n  x_1  geq 0\nendalign*One can verify that for any 0 leq x_2 leq 1, there exists a value x_1 such that (x_1 x_2) is in this polyhedron. This means that the H-representation obtained by eliminating x_1 is:beginalign*\n  x_1  leq 1 \n  x_1  geq 0\nendalign*where x_1 in the H-representation above represents x_2 in the previous one. This can be obtained as followsjulia> poly_x2 = eliminate(poly, [1])\njulia> hrep(poly_x2)\nH-representation\nbegin\n 2 2 rational\n 1//1 -1//1\n 0//1 1//1\nendThere is two methods of computing the elimination implemented in CDDLib: Fourier-Motzkin elimination and block elimination. As written by K. Fukuda in CDD's documentation, \"[Block elimination] might be a faster way to eliminate variables than the repeated [Fourier-Motzkin elimination] when the number of variables to eliminate is large\". You can specify the method to use as a third argument, e.g. eliminate(poly, [1], :FourierMotzkin), eliminate(poly, [1], :BlockElimination). A third method can be chosen: :ProjectGenerators. It computes the V-representation and then project each of its elements. This is the method of choice when the V-representation is already computed.If nothing is specified as in the block of code above, the behavior depends on the polyhedral library. If neither Fourier-Motzkin nor block elimination is implemented or if the V-representation is already computed then :ProjectGenerators is chosen. Otherwise, Polyhedra lets the library decide. In CDDLib, :FourierMotzkin is chosen when only the last dimension needs to be eliminated and :BlockElimination is chosen otherwise. Note that CDDLib only supports projecting the last trailing dimensions."
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
    "location": "utilities.html#Polyhedra.convexhull",
    "page": "Utilities",
    "title": "Polyhedra.convexhull",
    "category": "Function",
    "text": "convexhull(P1::VRep, P2::VRep)\n\nTakes the convex hull of P1 and P2  lambda x + (1-lambda) y  x in P1 y in P2 . It is very efficient between two V-representations or between two polyhedron for which the V-representation has already been computed. However, if P1 (resp. P2) is a polyhedron for which the V-representation has not been computed yet, it will trigger a representation conversion which is costly.\n\nThe type of the result will be chosen closer to the type of P1. For instance, if P1 is a polyhedron (resp. V-representation) and P2 is a V-representation (resp. polyhedron), convexhull(P1, P2) will be a polyhedron (resp. V-representation). If P1 and P2 are both polyhedra (resp. V-representation), the resulting polyhedron type (V-representation type) will be computing according to the type of P1. The coefficient type however, will be promoted as required taking both the coefficient type of P1 and P2 into account. ```\n\n\n\n"
},

{
    "location": "utilities.html#Operations-1",
    "page": "Utilities",
    "title": "Operations",
    "category": "section",
    "text": "convexhull"
},

{
    "location": "utilities.html#Polyhedra.ininterior",
    "page": "Utilities",
    "title": "Polyhedra.ininterior",
    "category": "Function",
    "text": "ininterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the interior of h. If h is an hyperplane, it always returns false. If h is an halfspace langle a x rangle leq beta, it returns whether p is in the open halfspace langle a x rangle  beta\n\nininterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the interior of h, e.g. in the interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.inrelativeinterior",
    "page": "Utilities",
    "title": "Polyhedra.inrelativeinterior",
    "category": "Function",
    "text": "inrelativeinterior(p::VRepElement, h::HRepElement)\n\nReturns whether p is in the relative interior of h. If h is an hyperplane, it is equivalent to p in h since the relative interior of an hyperplane is itself. If h is an halfspace, it is equivalent to ininterior(p, h).\n\ninrelativeinterior(p::VRepElement, h::HRep)\n\nReturns whether p is in the relative interior of h, e.g. in the relative interior of all the hyperplanes and halfspaces supporting h.\n\n\n\n"
},

{
    "location": "utilities.html#Containment-1",
    "page": "Utilities",
    "title": "Containment",
    "category": "section",
    "text": "ininterior\ninrelativeinterior"
},

{
    "location": "utilities.html#Polyhedra.volume",
    "page": "Utilities",
    "title": "Polyhedra.volume",
    "category": "Function",
    "text": "volume(p::Polyhedron{N, T}) where {N, T}\n\nReturns the N-dimensional hyper-volume of the polyhedron p. Returns Inf or -one(T) if it is infinite depending on whether the type T has an infinite value.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.surface",
    "page": "Utilities",
    "title": "Polyhedra.surface",
    "category": "Function",
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
    "category": "Function",
    "text": "chebyshevcenter(p::Rep[, solver])\n\nIf p is a H-representation or is a polyhedron for which the H-representation has already been computed, calls hchebyshevcenter, otherwise, call vchebyshevcenter.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.hchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.hchebyshevcenter",
    "category": "Function",
    "text": "hchebyshevcenter(p::HRep[, solver])\n\nReturn a tuple with the center and radius of the largest euclidean ball contained in the polyhedron p. Throws an error if the polyhedron is empty or if the radius is infinite.\n\n\n\n"
},

{
    "location": "utilities.html#Polyhedra.vchebyshevcenter",
    "page": "Utilities",
    "title": "Polyhedra.vchebyshevcenter",
    "category": "Function",
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
