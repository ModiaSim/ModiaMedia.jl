#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

dict["MoistAir"] = ModiaMedia.MoistAir(fluidData["H2O"],
                                       fluidData["N2"],
                                       singleGasesData["H2O"],
                                       singleGasesData["Air"])

