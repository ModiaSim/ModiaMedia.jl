module Runtests

using ModiaMedia
#using  Test


#@testset "Test ModiaMedia" begin
   Medium = getMedium("ConstantPropertyLiquidWater")

   println("\n... Medium = ", Medium)

   state  = setState_pT(Medium, 1e5, 300.0)
   d      = density(Medium,state)
   h      = specificEnthalpy(Medium,state)

   println("\nd = ", d)
   println("h = ", h)
   ModiaMedia.standardPlot(Medium; figure=1)


   Medium2 = getMedium("SimpleAir")

   println("\n... Medium2 = ", Medium2)

   state2  = setState_pT(Medium2, 1e5, 300.0)
   d2      = density(Medium2,state2)
   h2      = specificEnthalpy(Medium2,state2)

   println("\nd = ", d2)
   println("h   = ", h2)
   ModiaMedia.standardPlot(Medium2; figure=2)



   Medium3 = getMedium("N2")

   println("\n... Medium3 = ", Medium3)

   state3 = setState_pT(Medium3, 1e5, 300.0)
   d3     = density(Medium3,state3)
   h3     = specificEnthalpy(Medium3,state3)

   println("\nd = ", d3)
   println("h = ", h3)
   ModiaMedia.standardPlot(Medium3; figure=3)

#end

end