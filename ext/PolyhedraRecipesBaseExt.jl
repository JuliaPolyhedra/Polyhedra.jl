module PolyhedraRecipesBaseExt

import RecipesBase
import Polyhedra

RecipesBase.@recipe function f(p::Polyhedra.Polyhedron)
    seriestype --> :shape

    # Only hide the legend if the user DID NOT provide a custom label
    if !haskey(plotattributes, :label)
        legend --> false
    end

    # Explicitly force series_annotations to attach to the shape
    if haskey(plotattributes, :series_annotations)
        series_annotations := plotattributes[:series_annotations]
    end
   
    Polyhedra.planar_contour(p)
end

end
