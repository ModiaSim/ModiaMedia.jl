"""
    package ModiaMedia

Media models for use with Modia and other Julia packages.
It is planned that this package contains media property models to be used in package
[Modia](https://github.com/ModiaSim/Modia.jl), but also in other Julia packages.
The initial goal is to achieve a similar functionality
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media),
the standard media library for Modelica models, but with improvements based on Julia's features
such as multiple dispatch.

This package is currently under development.
"""
module ModiaMedia

const path    = dirname(dirname(@__FILE__))          # Absolute path of package directory
const Version = "0.1.0-dev from 2019-01-31 12:28"

println(" \nImporting ModiaMedia version ", Version)


export AbstractMedium, PureSubstance, getMedium
export MoistAir, SimpleMedium, SimpleIdealGasMedium, SingleGasNasa

export density, density_phX, density_pTX, density_der_1, density_pT, density_pT_der_1, density_pT_der_2, density_pT_der_3
export specificInternalEnergy_T, specificInternalEnergy_T_der_1, specificInternalEnergy_T_der_2

export temperature, temperature_phX, temperature_ph
export pressure, pressure_dT
export specificEnthalpy, specificEnthalpy_pTX, specificEnthalpy_dT, specificEnthalpy_T
export specificInternalEnergy, specificHeatCapacityCp
export setState_pTX, setState_pT, setState_ph, setState_ps, setState_dT
export isenthalpicState
export dynamicViscosity

# PureSubstance functions
export density_ph, temperature_ph, pressure_dT, specificEnthalpy_dT

# ThermodynamicState functions
export ThermodynamicState, ThermodynamicStates, ThermodynamicState_pT
export IndependentVariables, IndependentVariables_T, IndependentVariables_pT, IndependentVariables_ph
export IndependentVariables_phX, IndependentVariables_pT, IndependentVariables_dTX


### Abstract types -------------------------------------------------------------------------------
# The following structures and names are identical to Modelica.Media.Interfaces
# with the only exception, that "Partial" is replaced by "Abstract"

"`abstract type AbstractMedium` - Abstract type of all media"
abstract type AbstractMedium end

"`abstract type PureSubstance <: AbstractMedium` - Abstract type of all media consisting of a pure substance"
abstract type PureSubstance <: AbstractMedium end

"`abstract type MixtureMedium <: AbstractMedium` - Abstract type of all media consisting of a mixture"
abstract type MixtureMedium <: AbstractMedium end

"`abstract type CondensingGases <: AbstractMedium` - Abstract type of all media consisting of condensing media"
abstract type CondensingGases <: MixtureMedium end

"`abstract type ThermodynamicState` - Abstract type of all media states"
abstract type ThermodynamicState end

"`abstract type MixtureThermodynamicState <: ThermodynamicState` - Abstract type of the states of all media consisting of a mixture"
abstract type MixtureThermodynamicState <: ThermodynamicState end

"`abstract type AbstractFluidConstants` - Abstract type of all FluidConstants structures"
abstract type AbstractFluidConstants end


### Importing packages -------------------------------------------------------------------------------
using  JSON
using  StaticArrays
using  Unitful
import ModiaMath
import Serialization


### Including files for the ModiaMedia module --------------------------------------------------------
include("Interfaces/PartialMedium.jl")
include("Interfaces/PartialPureSubstance.jl")
include("Interfaces/PartialMixtureMedium.jl")
include("Interfaces/PartialCondensingGases.jl")

include("Media/SimpleMedium.jl")
include("Media/SimpleIdealGasMedium.jl")
include("Media/SingleGasNasa.jl")
include("Media/MoistAir.jl")


### Load medium dictionary from file
function loadMediumDict(file::AbstractString)
    mediumDict = Dict{AbstractMedium,Any}()
    if isfile(file)
        println("... Read media dictionary from file:\n",
                "    ", file)

        try
            f = open(file)
            mediumDict = Serialization.deserialize(f)
            close(f)
        catch
            error("\n\nFile \"", file, "\" is not compatible to the modified ModiaMedia source.\n",
                  "You have to delete this file and run `include(\"\$(ModiaMedia.path)/dict/GenerateMediumDict.jl\")` to regenerate it.\n")
        end

    else
        println("\n... File ", file, " does not exist.\n",
                "    You have to run `include(\"\$(ModiaMedia.path)/dict/GenerateMediumDict.jl\")` to generate it.\n",
                "    Without this file, no Medium can be used, because the Medium data is missing.")
    end
    return mediumDict
end

const mediumDictFile = "$path/src/Media/media.julia_serializer"
const mediumDict     = loadMediumDict(mediumDictFile)



### Inquire medium
getMedium(name::AbstractString) = length(mediumDict) == 0 ? loadMediumDict(mediumDictFile) : mediumDict[name]

end # module
