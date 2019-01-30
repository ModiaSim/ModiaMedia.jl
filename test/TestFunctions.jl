module TestFunctions

using ModiaMedia
using Test

const mediaNames = ["ConstantPropertyLiquidWater", "SimpleAir", "N2"]

function testMediumFunctions(mediumName; figure=1)
   Medium = getMedium(mediumName)
   infos  = Medium.infos
   state  = setState_pTX(Medium, infos.reference_p,
                                 infos.reference_T, 
                                 infos.reference_X)

   p   = pressure(state)
   T   = temperature(state)
   h   = specificEnthalpy(state)
   d   = density(state)
   u   = specificInternalEnergy(state)
   Cp  = specificHeatCapacityCp(state)
   eta = dynamicViscosity(state)
 
   println("... ", mediumName, " (figure=",figure,"): p = ", p, ", T = ", T, ", h = ", h, ", d = ", d, ", u = ", u, ", Cp = ", Cp, ", eta = ", eta)

   state_b = isenthalpicState(state, 0.1e5)
   h_b     = specificEnthalpy(state_b)
   @test isapprox(h,h_b)

   if typeof(Medium) <: PureSubstance
      d2 = density_ph(Medium,p,h)
      T2 = temperature_ph(Medium,p,h)
      @test isapprox(d,d2)
      @test isapprox(T,T2)

      if !Medium.infos.singleState
         p2 = pressure_dT(Medium,d,T)
         h2 = specificEnthalpy_dT(Medium,d,T)
         @test isapprox(p,p2)
         @test isapprox(h,h2)
      end
   end

   ModiaMedia.standardPlot(Medium; figure=figure)
end


function testAllMediumFunctions(mediaNames)
   println("\n... Test thermodynamic property functions:")
   i = 0
   for name in mediaNames
      i += 1
      testMediumFunctions(name; figure=i)
   end
end

@testset "ModiaMedia/test/TestFunctions.jl:" begin 
   testAllMediumFunctions(mediaNames)
end

end
