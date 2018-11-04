#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

const simpleMediaDict = JSON.parsefile(joinpath(dirname(@__FILE__), "SimpleMedia.json"))


function SimpleMedium(name::AbstractString, data)
    medium = simpleMediaDict[name]    

    ModiaMedia.SimpleMedium(infos     = fillobj(medium["infos"]    , ModiaMedia.FluidInfos(mediumName=name, singleState=true,
                                                                                           ThermoStates=ModiaMedia.IndependentVariables_T) ),
                            constants = fillobj(medium["constants"], ModiaMedia.BasicFluidConstants() ),
                            limits    = fillobj(medium["limits"]   , ModiaMedia.FluidLimits() ),
                            data      = fillobj(medium["data"]     , ModiaMedia.SimpleMediumData() )
                            )
end


function storeSimpleMedium!(mediumDict)
    for (name,data) in simpleMediaDict
        m = SimpleMedium(name, data)
        p = m.infos.p_default
        m.infos.h_default = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.infos.T_default))
        m.limits.HMIN     = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.limits.TMIN))
        m.limits.HMAX     = ModiaMedia.specificEnthalpy(m, ModiaMedia.setState_pT(m,p,m.limits.TMAX))
        mediumDict[name]  = m
    end
end

storeSimpleMedium!(dict)
