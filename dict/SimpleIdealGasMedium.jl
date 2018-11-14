#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

dict["SimpleAir"] = ModiaMedia.SimpleIdealGasMedium(
                        mediumName = "SimpleAir",
                        fluidConstants = ModiaMedia.BasicFluidConstants(
                            iupacName="simple air",
                            casRegistryNumber="not a real substance",
                            chemicalFormula="N2, O2",
                            structureFormula="N2, O2",
                            molarMass=0.0289651159),
                        data = ModiaMedia.SimpleIdealGasMediumData(
                            cp_const     = 1005.45, 
                            MM_const     = 0.0289651159, 
                            R_gas        = 8.3144598/0.0289651159, 
                            eta_const    = 1.82e-5, 
                            lambda_const = 0.026, 
                            T_min        = ModiaMedia.from_degC(0.0), 
                            T_max        = ModiaMedia.from_degC(100.0),
                            T0           = 298.15)
                    )

