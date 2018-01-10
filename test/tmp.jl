using Polyhedra, CDDLib
P = polyhedron(SimpleVRepresentation(eye(2)), CDDLibrary())
P * P
