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
  # (could be also set from p,h, or p,s, or d,T, or
  # p,T,X, or p,h,X, or p,s,X, or d,T,X)
  state = setState_pT(Medium, p, T)

  # Update a state object with new values
  setState_pT!(state, p, T)

  # Call media functions (here to compute density and specific enthalpy)
  d = density(state)
  h = specificEnthalpy(state)

  # Print computed values
  println("data for p=$p, T=$T:")
  println("density          = ", d)
  println("specificEnthalpy = ", h)

  # List the available media
  listMedia()

  # Plot the most important characteristics of the medium
  standardPlot(Medium)
```

This example generates the following plot:

![standardPlot](../resources/images/N2.png)



## Available media

The available media can be listed with `listMedia()` resulting in:

| Row | name                        | type                 |
|-----|-----------------------------|----------------------|
| 1   | Ar                          | [SingleGasNasa](@ref)        |
| 2   | C2H2\_vinylidene             | [SingleGasNasa](@ref)        |
| 3   | C2H4                        | [SingleGasNasa](@ref)        |
| 4   | C2H5OH                      | [SingleGasNasa](@ref)        |
| 5   | C2H6                        | [SingleGasNasa](@ref)        |
| 6   | C3H6\_propylene              | [SingleGasNasa](@ref)        |
| 7   | C3H8                        | [SingleGasNasa](@ref)        |
| 8   | C4H10\_n\_butane              | [SingleGasNasa](@ref)        |
| 9   | C4H8\_1\_butene               | [SingleGasNasa](@ref)        |
| 10  | C5H10\_1\_pentene             | [SingleGasNasa](@ref)        |
| 11  | C5H12\_n\_pentane             | [SingleGasNasa](@ref)        |
| 12  | C6H12\_1\_hexene              | [SingleGasNasa](@ref)        |
| 13  | C6H14\_n\_hexane              | [SingleGasNasa](@ref)        |
| 14  | C6H6                        | [SingleGasNasa](@ref)        |
| 15  | C7H14\_1\_heptene             | [SingleGasNasa](@ref)        |
| 16  | C7H16\_n\_heptane             | [SingleGasNasa](@ref)        |
| 17  | C8H10\_ethylbenz             | [SingleGasNasa](@ref)        |
| 18  | C8H18\_n\_octane              | [SingleGasNasa](@ref)        |
| 19  | CH3OH                       | [SingleGasNasa](@ref)        |
| 20  | CH4                         | [SingleGasNasa](@ref)        |
| 21  | CL2                         | [SingleGasNasa](@ref)        |
| 22  | CO                          | [SingleGasNasa](@ref)        |
| 23  | CO2                         | [SingleGasNasa](@ref)        |
| 24  | ConstantPropertyLiquidWater | [SimpleMedium](@ref)         |
| 25  | F2                          | [SingleGasNasa](@ref)        |
| 26  | H2                          | [SingleGasNasa](@ref)        |
| 27  | H2O                         | [SingleGasNasa](@ref)        |
| 28  | He                          | [SingleGasNasa](@ref)        |
| 29  | MoistAir                    | [MoistAir](@ref)             |
| 30  | N2                          | [SingleGasNasa](@ref)        |
| 31  | N2O                         | [SingleGasNasa](@ref)        |
| 32  | NH3                         | [SingleGasNasa](@ref)        |
| 33  | NO                          | [SingleGasNasa](@ref)        |
| 34  | NO2                         | [SingleGasNasa](@ref)        |
| 35  | Ne                          | [SingleGasNasa](@ref)        |
| 36  | O2                          | [SingleGasNasa](@ref)        |
| 37  | SO2                         | [SingleGasNasa](@ref)        |
| 38  | SO3                         | [SingleGasNasa](@ref)        |
| 39  | SimpleAir                   | [SimpleIdealGasMedium](@ref) |


## Structure of package

A medium is a struct of the following type:

```julia
mutable struct MediumName <: ModiaMedia.AbstractMedium  # or of a subtype of AbstractMedium
    infos::ModiaMedia.FluidInfos
    fluidConstants::Vector{ModiaMedia.AbstractFluidConstants}
    fluidLimits::ModiaMedia.FluidLimits
    data  # medium specific data
end

struct FluidInfos
    mediumName::AbstractString                   # "Name of the medium";
    substanceNames::Vector{AbstractString}       # "Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.";
    extraPropertiesNames::Vector{AbstractString} # "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused"
    ThermoStates::IndependentVariables           # "Enumeration type for independent variables";
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



## Main Developers

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/) ([DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en))\
- Hilding Elmqvist ([Mogram](http://www.mogram.net/)),\
- [Chris Laughman](http://www.merl.com/people/laughman) ([MERL](http://www.merl.com/)).
- All the content of ModiaMedia is based on
  [Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media)
  which was and is developed from many people.

License: MIT (expat)


## Release Notes

The ModiaMedia package development has just started and a lot has to be improved.


### Version 0.1.0-dev

A version is not yet released.
