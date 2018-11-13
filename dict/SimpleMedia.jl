#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

const simpleMediaDict = JSON.parsefile(joinpath(dirname(@__FILE__), "SimpleMedia.json"))

function SimpleMedium(name::AbstractString, data)
    medium = simpleMediaDict[name]    

    ModiaMedia.SimpleMedium(mediumName  = medium["mediumName"],
                            reference_p = medium["reference_p"],
                            reference_T = medium["reference_T"],
                            p_default   = medium["p_default"],
                            T_default   = medium["T_default"],
                            fluidConstants = fillobj(medium["fluidConstants"], ModiaMedia.BasicFluidConstants() ),
                            fluidLimits    = fillobj(medium["fluidLimits"]   , ModiaMedia.FluidLimits() ),
                            data           = fillobj(medium["data"]          , ModiaMedia.SimpleMediumData() )
                            )
end


function storeSimpleMedium!(mediumDict)
    for (name,value) in simpleMediaDict
        mediumDict[name] = SimpleMedium(name, value)
    end
end

storeSimpleMedium!(dict)
