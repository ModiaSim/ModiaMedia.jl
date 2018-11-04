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
const Version = "0.1.0-dev from 2018-11-04 18:34"

println(" \nImporting ModiaMedia version ", Version)



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
include("Media/SingleGasNasa.jl")


### Load medium dictionary from file
const file = "$path/src/Media/media.julia_serializer"
println("... Read media dictionary from file:\n",
        "    ", file)

f = open(file)  
const mediumDict = Serialization.deserialize(f)
close(f)
 
### Inquire medium
Medium(name::AbstractString) = mediumDict[name]


end # module
