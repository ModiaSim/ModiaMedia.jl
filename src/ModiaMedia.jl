"""
    package ModiaMedia

Media models for use with Modia and other Julia packages.
It is planned that this package contains media property models to be used in package
[Modia](https://github.com/ModiaSim/Modia.jl), but also in other Julia packages.
The initial goal is to achieve a similar functionality as
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media),
the standard media library for Modelica models, but with improvements based on Julia's features
such as multiple dispatch.

A medium model is basically a struct that has the following structure:

```
mutable struct MediumName <: AbstractMedium
    infos::ModiaMedium.FluidInfos
    fluidConstants::SVector{1,ModiaMedium.AbstractFluidConstants}
    fluidLimits::ModiaMedium.FluidLimits
    data::MediumSpecificData
end
```

and functions operating on instances of such a struct.

This package is currently under development.
"""
module ModiaMedia

const path    = dirname(dirname(@__FILE__))          # Absolute path of package directory
const Version = "0.1.0-dev"
const Date    = "2019-02-25"

println(" \nImporting ModiaMedia Version $Version ($Date)")



# Export Abstract types
export AbstractMedium
export PureSubstance
export MixtureMedium
export CondensingGases
export ThermodynamicState
export MixtureThermodynamicState
export AbstractFluidConstants

# Export types used in a Medium
export FluidInfos
export FluidLimits
export BasicFluidConstants
export IdealGasFluidConstants


# Export State structs
export SimpleMediumState
export SingleGasNasaState
export SimpleIdealGasMediumState
export MoistAirState


# Export general functions
export getMedium, listMedia, standardCharacteristics, standardPlot


# Export set state functions
export setState_pTX, setState_pTX!, setState_pT, setState_pT!
export setState_phX, setState_phX!, setState_ph, setState_ph!
export setState_psX, setState_psX!, setState_ps, setState_ps!
export setState_dTX, setState_dTX!, setState_dT, setState_dT!
export isenthalpicState, isenthalpicState!


# Export thermodynamic property functions for all media
export density, density_phX, density_pTX, density_der_1, density_pT, density_pT_der_1, density_pT_der_2, density_pT_der_3
export specificInternalEnergy_T, specificInternalEnergy_T_der_1, specificInternalEnergy_T_der_2
export temperature, temperature_phX, temperature_ph
export pressure, pressure_dT
export specificEnthalpy, specificEnthalpy_pTX, specificEnthalpy_dT, specificEnthalpy_T
export specificInternalEnergy, specificHeatCapacityCp
export dynamicViscosity
export gasConstant


# Export PureSubstance thermodynamic property  functions
export density_ph, temperature_ph, pressure_dT, specificEnthalpy_dT


# Export MoistAir specific functions
export saturationPressureOfLiquidWater
export sublimationPressureIce
export saturationPressureOfWater


# Enumeration values
export IndependentVariables, IndependentVariables_T, IndependentVariables_pT, IndependentVariables_ph
export IndependentVariables_phX, IndependentVariables_pT, IndependentVariables_dTX
export ReferenceEnthalpy, ReferenceEnthalpy_ZeroAt0K, ReferenceEnthalpy_ZeroAt25C, ReferenceEnthalpy_UserDefined
export ReferenceEntropy , ReferenceEntropy_ZeroAt0K , ReferenceEntropy_ZeroAt0C  , ReferenceEntropy_UserDefined



### Abstract types -------------------------------------------------------------------------------
# The following structures and names are identical to Modelica.Media.Interfaces
# with the only exception, that "Partial" is replaced by "Abstract"

"""
    abstract type AbstractMedium

Abstract type of all media.
"""
abstract type AbstractMedium end


"""
    abstract type PureSubstance <: AbstractMedium

Abstract type of all media consisting of a pure substance.
"""
abstract type PureSubstance <: AbstractMedium end


"""
    abstract type MixtureMedium <: AbstractMedium

Abstract type of all media consisting of a mixture of media.
"""
abstract type MixtureMedium <: AbstractMedium end


"""
    abstract type CondensingGases <: AbstractMedium

Abstract type of all media consisting of condensing media.
"""
abstract type CondensingGases <: MixtureMedium end


"""
   abstract type ThermodynamicState

Abstract type of all media states.
"""
abstract type ThermodynamicState end


"""
    abstract type MixtureThermodynamicState <: ThermodynamicState

Abstract type of the states of all media consisting of a mixture of media.
"""
abstract type MixtureThermodynamicState <: ThermodynamicState end


"""
    abstract type AbstractFluidConstants

Abstract type of all FluidConstants structures.
"""
abstract type AbstractFluidConstants end



### Importing packages -------------------------------------------------------------------------------
using  JSON
using  StaticArrays
using  Unitful
import DataFrames
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

const mediumDictFile   = "$path/src/Media/media.julia_serializer"
const generateDictFile = "$path/dict/GenerateMediumDict.jl"
const mediumDict       = Array{Any}(nothing,1)


### Load medium dictionary from file
function loadMediumDictFile()
    global mediumDictFile
    println("... Read media dictionary from file:\n",
            "    ", mediumDictFile)
    f = open(mediumDictFile)
    dict = Serialization.deserialize(f)
    close(f)
    return dict
end

function loadMediumDict()
    global generateDictFile

    if isfile(mediumDictFile)
        try
            dict = loadMediumDictFile()
        catch
            println("... Regenerate media dictionary")
            include(generateDictFile)
            dict = loadMediumDictFile()
        end
    else
        include(generateDictFile)
        dict = loadMediumDictFile()
    end

    return dict
end



"""
    Medium = getMedium(name::AbstractString)

Return `Medium` object from medium `name`. 
Possible values of argument `name` can be inquired via `listMedia()`
([Available media](@ref)).

# Examples
```julia
medium = getMedium("SimpleAir")
```

"""
function getMedium(name::AbstractString)::AbstractMedium
    global mediumDict
    if typeof(mediumDict[1]) == Nothing
        mediumDict[1] = loadMediumDict()
    end
    return mediumDict[1][name]
end


"""
    listMedia()

List available media of ModiaMedia.
"""
function listMedia()::Nothing
    global mediumDict
    if typeof(mediumDict[1]) == Nothing
        mediumDict[1] = loadMediumDict()
    end

    dict = mediumDict[1]
    media_table = DataFrames.DataFrame(name=AbstractString[], type=Symbol[])

    for key in sort(collect(keys(dict)))
        push!(media_table, [key, Symbol(typeof(dict[key]))])
    end

    println("\nMedia available in ModiaMedia:\n")
    show(media_table,allrows=true,allcols=true,summary=false)
    return nothing
end



# Import packages that are used in examples and tests
# (in order that there are no requirements on the environment
#  in which the examples and tests are executed).
import JSON
import ModiaMath
import Serialization
import StaticArrays
import Unitful
import Test




end # module