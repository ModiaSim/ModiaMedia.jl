#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

const simpleMediaDict = JSON.parsefile(joinpath(dirname(@__FILE__), "SimpleMedia.json"))

function SimpleMedium(name::AbstractString, data)
    medium = simpleMediaDict[name]    

    ModiaMedia.SimpleMedium(infos          = fillobj(medium["infos"]         , ModiaMedia.FluidInfos(mediumName=name, singleState=true,
                                                                                  baseProperties= :BaseProperties_SimpleMedium,
                                                                                  ThermoStates=ModiaMedia.IndependentVariables_T) ),
                            fluidConstants = fillobj(medium["fluidConstants"], ModiaMedia.BasicFluidConstants() ),
                            fluidLimits    = fillobj(medium["fluidLimits"]   , ModiaMedia.FluidLimits() ),
                            data           = fillobj(medium["data"]          , ModiaMedia.SimpleMediumData() )
                            )
end


function storeSimpleMedium!(mediumDict)
    for (name,data) in simpleMediaDict
        m = SimpleMedium(name, data)
        p = m.infos.p_default
        m.infos.h_default = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.infos.T_default))
        m.fluidLimits.HMIN = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.fluidLimits.TMIN))
        m.fluidLimits.HMAX = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.fluidLimits.TMAX))
        mediumDict[name]    = m
    end
end

storeSimpleMedium!(dict)
