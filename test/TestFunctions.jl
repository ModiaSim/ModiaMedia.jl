module TestFunctions

using ModiaMedia
using ModiaMedia.Test  # included via ModiaMedia, to avoid requirement to add it in the standard environment  


const mediaNames = ["ConstantPropertyLiquidWater", "SimpleAir", "N2", "MoistAir"]

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
   #Cp  = specificHeatCapacityCp(state)
   #eta = dynamicViscosity(state)
 
   #println("... ", mediumName, " (figure=",figure,"): p = ", p, ", T = ", T, ", h = ", h, ", d = ", d, ", u = ", u, ", Cp = ", Cp, ", eta = ", eta)
   println("... ", mediumName, " (figure=",figure,"): p = ", p, ", T = ", T, ", h = ", h, ", d = ", d, ", u = ", u)

   state_b = isenthalpicState(state, 0.1e5)
   h_b     = specificEnthalpy(state_b)
   @test isapprox(h,h_b)

   state_b = setState_phX(Medium, p, h, infos.reference_X)
   T_b     = temperature(state_b)
   @test isapprox(T,T_b)
   
   if !infos.singleState
      state_b = setState_dTX(Medium, d, T, infos.reference_X)
      h_b     = specificEnthalpy(state_b)
      @test isapprox(h,h_b)
   end

   state_c = setState_pTX(Medium, 0.9*infos.reference_p,
                                  0.9*infos.reference_T, 
                                  0.9*infos.reference_X)
   state_d = deepcopy(state_c)
   setState_pTX!(state_d, p,T,infos.reference_X)
   @test isapprox(p, pressure(state_d))
   @test isapprox(T, temperature(state_d))

   state_d = deepcopy(state_c)
   setState_phX!(state_d, p,h,infos.reference_X)
   @test isapprox(p, pressure(state_d))
   @test isapprox(h, specificEnthalpy(state_d))

   if !infos.singleState
      state_d = deepcopy(state_c)
      setState_dTX!(state_d, d,T,infos.reference_X)
      @test isapprox(d, density(state_d))
      @test isapprox(T, temperature(state_d))
   end

   if typeof(Medium) <: PureSubstance
      state_b = setState_pT(Medium, p, T)
      T_b     = temperature(state_b)
      @test isapprox(T,T_b)

      state_b = setState_ph(Medium, p, h)
      T_b     = temperature(state_b)
      @test isapprox(T,T_b)

      state_d = deepcopy(state_c)
      setState_pT!(state_d, p,T)
      @test isapprox(p, pressure(state_d))
      @test isapprox(T, temperature(state_d))

      state_d = deepcopy(state_c)
      setState_ph!(state_d, p,h)
      @test isapprox(p, pressure(state_d))
      @test isapprox(h, specificEnthalpy(state_d))

      d2 = density_ph(Medium,p,h)
      T2 = temperature_ph(Medium,p,h)
      @test isapprox(d,d2)
      @test isapprox(T,T2)

      if !infos.singleState
         state_b = setState_dT(Medium, d, T)
         h_b     = specificEnthalpy(state_b)
         @test isapprox(h,h_b)

         state_d = deepcopy(state_c)
         setState_dT!(state_d, d,T)
         @test isapprox(d, density(state_d))
         @test isapprox(T, temperature(state_d))

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
