module Runtests

import ModiaMedia
#using  Test


#@testset "Test ModiaMedia" begin
   medium = ModiaMedia.Medium("SimpleLiquidWater")

   println("\n... medium = ", medium)

   state  = ModiaMedia.setState_pT(medium, 1e5, 300.0)
   d      = ModiaMedia.density(medium,state)
   h      = ModiaMedia.specificEnthalpy(medium,state)

   println("\nd = ", d)
   println("h = ", h)
   ModiaMedia.standardPlot(medium; figure=1)


   medium2 = ModiaMedia.Medium("N2")

   println("\n... medium2 = ", medium2)

   state2 = ModiaMedia.setState_pT(medium2, 1e5, 300.0)
   d2     = ModiaMedia.density(medium2,state2)
   h2     = ModiaMedia.specificEnthalpy(medium2,state2)

   println("\nd = ", d2)
   println("h = ", h2)
   ModiaMedia.standardPlot(medium2; figure=2)

#end

end