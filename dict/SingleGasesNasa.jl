#
# This file is part of module 
#   GenerateMediumDict (ModiaMedia/dict/GenerateMediumDict.jl)
#

const Modelica_Media_IdealGases_Common_DataRecord = ModiaMedia.SingleGasNasaData
const singleGasesData = Dict{AbstractString, ModiaMedia.SingleGasNasaData}()
const fluidData       = Dict{AbstractString, ModiaMedia.IdealGasFluidConstants}()


# Include data from Modelica.Media.IdealGases.Common.SingleGasesData.jl into dict
include(joinpath(dirname(@__FILE__), "Modelica.Media.IdealGases.Common.SingleGasesData.jl"))

# Include data from Modelica.Media.IdealGases.Common.FluidData.jl into dict
include(joinpath(dirname(@__FILE__), "Modelica.Media.IdealGases.Common.FluidData.jl"))


# Construct Medium
"""
    medium = SimpleGasNasa(name, fluidConstants, data)

Return an ideal gas medium of a single substance, or an ideal gas medium 
consisting of several substances where the mass fractions are fixed. 
Density is a function of T and p. All other quantities are solely a function of T.


# Sources for model and literature

Original Data: *Computer program for calculation of complex chemical equilibrium compositions and applications*.
Part 1: Analysis Document ID: 19950013764 N (95N20180) 
File Series: NASA Technical Reports Report 
Number: NASA-RP-1311 E-8017 NAS 1.61:1311 
Authors: Gordon, Sanford (NASA Lewis Research Center) Mcbride, Bonnie J. (NASA Lewis Research Center) 
Published: Oct 01, 1994. 

The data to computpe the specific enthalpy is from
McBride B.J., Zehe M.J., and Gordon S. (2002): 
*NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species*. 
NASA report TP-2002-211556.
https://www.grc.nasa.gov/WWW/CEAWeb/TP-2002-211556.pdf

Known limits of validity: The data is valid for temperatures between 200K and 6000K.
A few of the data sets for monatomic gases have a discontinuous 1st derivative at 1000K, 
but this never caused problems so far. 

This model has been adapted from the Modelica.Media.IdealGases.Common.SingleGasNasa 
Modelica package, which in turn was copied from the ThermoFluid library and 
adapted to Modelica. 
"""
function SingleGasNasa(name::AbstractString, fluidConstants::ModiaMedia.IdealGasFluidConstants, data::ModiaMedia.SingleGasNasaData)
    ModiaMedia.SingleGasNasa(mediumName     = name,
                             fluidConstants = fluidConstants,
                             fluidLimits    = ModiaMedia.FluidLimits(TMIN=200.0, TMAX=6000.0),
                             data           = data)
end


# Construct SingleGasNasa objects and store them in the medium dict
function storeSingleGasNasaMedium!(mediumDict)
    global fluidData
    global singleGasesData
    for (name,fluidConstants) in fluidData
        dict[name] = SingleGasNasa(name, fluidConstants, singleGasesData[name])
    end
end

# Manually insert air into the dictionary according to model from Modelica.Media.Air.MoistAir; should probably fix somehow 
dict["Air"] = SingleGasNasa("Air", fluidData["N2"], singleGasesData["Air"])

storeSingleGasNasaMedium!(dict)
