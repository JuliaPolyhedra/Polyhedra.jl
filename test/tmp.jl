using Polyhedra, CDDLib
P = polyhedron(vrep(eye(2)), CDDLibrary())
P * P
