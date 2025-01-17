module PolyhedraRecipesBaseExt

import RecipesBase
import Polyhedra

RecipesBase.@recipe function f(p::Polyhedra.Polyhedron)
    seriestype --> :shape
    legend --> false
    Polyhedra.planar_contour(p)
end

end
