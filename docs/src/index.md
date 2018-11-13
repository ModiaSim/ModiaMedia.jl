# ModiaMedia.jl Documentation

[ModiaMedia](https://github.com/ModiaSim/ModiaMedia.jl) shall provide Media models 
for use with [Modia](https://github.com/ModiaSim/Modia.jl)
and other Julia packages. The initial goal is to achieve a similar functionality as
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media),
the standard media library for Modelica models, but with improvements based on Julia features
such as multiple dispatch.

This package is under development and it is planned to provide all media from
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media)
in this package.


## Installation

This package is currently under development and is not yet registered in METADATA.
Julia 1.0 is required. Installation is performed via:

```julia
julia> ]add https://github.com/ModiaSim/ModiaMedia.jl
```

ModiaMedia uses [PyPlot](https://github.com/JuliaPy/PyPlot.jl) for plotting (via ModiaMath.plot).
If `PyPlot` is not available in your current Julia environment
an information message is printed and all `plot(..)` calls are ignored.

In order that plot windows are displayed, you need to add `PyPlot` to your current environment
via `]add PyPlot`. Often this automatic installation fails and it is recommended to follow
the instructions
[Installing PyPlot in a robust way](https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way).


## Use

```julia
  using ModiaMedia

  # Define medium to be used
  Medium = getMedium("N2");

  # Define the operating point where the medium shall be evaluated.
  p = 1e5    # in [Pa]
  T = 300.0  # in [K]

  # Set the medium-specific thermodynamic state from p and T
  # (could be also set from p and h, or p and s, or d and T)
  state = setState_pT(Medium, p, T)

  # Call media functions (here to compute density and specific enthalpy)
  d = density(Medium,state)
  h = specificEnthalpy(Medium,state)

  # Print computed values
  println("data for p=$p, T=$T:")
  println("density          = ", d)
  println("specificEnthalpy = ", h)

  # Plot the most important characteristics of the medium
  ModiaMedia.standardPlot(Medium)
```

This example generates the following plot:

![standardPlot](../resources/images/N2.png)



### Currently available media

- SimpleLiquidWater

- The following 37 ideal gases (from NASA Glenn coefficients): 
  Ar, CH4, CH3OH, CO, CO2, C2H2_vinylidene, C2H4, C2H5OH, C2H6, C3H6_propylene, 
  C3H8, C3H8O_1propanol, C4H8_1_butene, C4H10_n_butane, C5H10_1_pentene, 
  C5H12_n_pentane, C6H6, C6H12_1_hexene, C6H14_n_hexane, C7H14_1_heptene,  
  C7H16_n_heptane, C8H10_ethylbenz, C8H18_n_octane, CL2, F2, H2, H2O,    
  He, NH3, NO, NO2, N2, N2O, Ne, O2, SO2, SO3 


### Structure of package

A medium is a struct of the following type:

```julia
struct MediumXXX <: AbstractMedium  # or of a subtype of AbstractMedium
    infos::FluidInfos
    fluidConstants::Vector{AbstractFluidConstants}
    fluidLimits::FluidLimits
    data  # medium specific data
end

struct FluidInfos 
    mediumName::AbstractString                   # "Name of the medium";
    substanceNames::Vector{AbstractString}       # "Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.";
    extraPropertiesNames::Vector{AbstractString} # "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused"
    ThermoStates::IndependentVariables           # "Enumeration type for independent variables";
    baseProperties::Symbol                       # "Symbol of baseProperties model = :BaseProperties_<StructName>
    singleState::Bool                            # "= true, if u and d are not a function of pressure";
    reducedX::Bool                               # "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
    fixedX::Bool                                 # "= true if medium contains the equation X = reference_X";
    reference_p::Float64                         # "Reference pressure of Medium: default 1 atmosphere";
    reference_T::Float64                         # "Reference temperature of Medium: default 25 deg Celsius";
    reference_X::AbstractVector                  # "Default mass fractions of medium";
    p_default::Float64                           # "Default value for pressure of medium (for initialization)";
    T_default::Float64                           # "Default value for temperature of medium (for initialization)";
    h_default::Float64                           # "Default value for specific enthalpy of medium (for initialization)";
    X_default::Vector{Float64}                   # "Default value for specific enthalpy of medium (for initialization)";
    nS::Int                                      # "Number of substances"
    nX::Int                                      # "Number of mass fractions"
    nXi::Int                                     # "Default value for mass fractions of medium (for initialization)"
    nC::Int                                      # "Number of extra (outside of standard mass-balance) transported properties"
    C_nominal::Vector{Float64}                   # "Default for the nominal values for the extra properties"  
end

struct BasicFluidConstants <: AbstractFluidConstants
    iupacName::AbstractString           # "Complete IUPAC name (or common name, if non-existent)";
    casRegistryNumber::AbstractString   # "Chemical abstracts sequencing number (if it exists)";
    chemicalFormula::AbstractString     # "Chemical formula, (brutto, nomenclature according to Hill";
    structureFormula::AbstractString    # "Chemical structure formula";
    molarMass::Float64                  # "Molar mass";
end

mutable struct FluidLimits
    TMIN::Float64  # "Minimum temperature";
    TMAX::Float64  # "Maximum temperature";
    DMIN::Float64  # "Minimum density";
    DMAX::Float64  # "Maximum density";
    PMIN::Float64  # "Minimum pressure";
    PMAX::Float64  # "Maximum pressure";
    HMIN::Float64  # "Minimum enthalpy";
    HMAX::Float64  # "Maximum enthalpy";
    SMIN::Float64  # "Minimum entropy";
    SMAX::Float64  # "Maximum entropy";
end
```

and all instances of this struct are stored in a **dictionary**.
This dictionary is constructed in a preprocessing step
by running "ModiaMedia/dict/GenerateMediumDict.jl".
This module contains code that was mostly automatically
converted from Modelica.Media to Julia.
The resulting dictionary is serialized and stored in "ModiaMedia/src/Media/media.julia_serializer".
When package ModiaMedia is compiled, this serialized dictionary is deserialized
and included in the compiled package.

Function `ModiaMedia.Medium(name)` returns the `MediumXXX` instance stored
in the medium dictionary with key `name`.


### Status

The ModiaMedia package development has just started and a lot has to be improved.


## Release Notes

### Version 0.1.0-dev

A version is not yet released.
