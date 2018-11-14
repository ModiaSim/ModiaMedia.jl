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
const Version = "0.1.0-dev from 2018-11-14 11:00"

println(" \nImporting ModiaMedia version ", Version)


export getMedium
export density, pressure, specificEnthalpy, specificInternalEnergy, temperature
export setState_pT, setState_ph, setState_ps, setState_dT


### Abstract types -------------------------------------------------------------------------------
# The following structures and names are identical to Modelica.Media.Interfaces
# with the only exception, that "Partial" is replaced by "Abstract"

"`abstract type AbstractMedium` - Abstract type of all media"
abstract type AbstractMedium end

"`abstract type PureSubstance <: AbstractMedium` - Abstract type of all media consisting of a pure substance"
abstract type PureSubstance <: AbstractMedium end
 
"`abstract type ThermodynamicState` - Abstract type of all media states"
abstract type ThermodynamicState end

"`abstract type PureSubstanceThermodynamicState <: ThermodynamicState` - Abstract type of the states of all media consisting of a pure substance"
abstract type PureSubstanceThermodynamicState <: ThermodynamicState end

"`abstract type AbstractFluidConstants` - Abstract type of all FluidConstants structures"
abstract type AbstractFluidConstants end



### Importing packages -------------------------------------------------------------------------------
using  JSON
using  StaticArrays
using  Unitful
import ModiaMath
import Serialization

 
### Including files for the ModiaMedia module --------------------------------------------------------
include("Interfaces/Unitful_U_str.jl")   # only temporarily until the u".." issue is fixed in Unitful
using  .Unitful_U_str

include("Interfaces/PartialMedium.jl")
include("Interfaces/PartialPureSubstance.jl")
include("Media/SimpleMedium.jl")
include("Media/SimpleIdealGasMedium.jl")
include("Media/SingleGasNasa.jl")


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
getMedium(name::AbstractString) = mediumDict[name]


end # module
