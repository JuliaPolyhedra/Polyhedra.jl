using Polyhedra
using Base.Test
using CDDLib
using GeometryTypes

ext1 = SimpleVRepresentation([0 1 2 3;0 2 1 3; 1 0 2 3; 1 2 0 3; 2 0 1 3; 2 1 0 3;
                              0 1 3 2;0 3 1 2; 1 0 3 2; 1 3 0 2; 3 0 1 2; 3 1 0 2;
                              0 3 2 1;0 2 3 1; 3 0 2 1; 3 2 0 1; 2 0 3 1; 2 3 0 1;
                              3 1 2 0;3 2 1 0; 1 3 2 0; 1 2 3 0; 2 3 1 0; 2 1 3 0])
poly1 = CDDPolyhedron{4,Rational{BigInt}}(ext1)

poly2 = project(poly1, [1 1 1; -1 1 1; 0 -2 1; 0 0 -3])
