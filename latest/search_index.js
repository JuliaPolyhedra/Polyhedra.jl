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
    "text": "Pages = [\"installation.md\", \"representation.md\", \"polyhedron.md\"]\nDepth = 2"
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

]}
