#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

dict["ConstantPropertyLiquidWater"] = ModiaMedia.SimpleMedium(
                        mediumName = "ConstantPropertyLiquidWater",
                        fluidConstants = ModiaMedia.BasicFluidConstants(
                            chemicalFormula="H2O", 
                            structureFormula="H2O", 
                            casRegistryNumber="7732-18-5", 
                            iupacName="oxidane", 
                            molarMass=0.018015268),
                        data = ModiaMedia.SimpleMediumData(
                            cp_const=4184, 
                            cv_const=4184, 
                            d_const=995.586, 
                            eta_const=1.e-3, 
                            lambda_const=0.598, 
                            a_const=1484, 
                            T_min=ModiaMedia.from_degC(-1), 
                            T_max=ModiaMedia.from_degC(130), 
                            T0=273.15, 
                            MM_const=0.018015268)
                    )
