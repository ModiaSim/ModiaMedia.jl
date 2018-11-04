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
  import ModiaMedia

  # Define medium to be used
  medium = ModiaMedia.Medium("N2");

  # Define the operating point where the medium shall be evaluated. 
  p = 1e5    # in [Pa]
  T = 300.0  # in [K]

  # Determine the medium-specific thermodynamic state
  # (setState_ph, setState_ps, setState_dT could also be used)
  state = ModiaMedia.setState_pT(medium, p, T)

  # Call media functions (here to compute density and specific enthalpy)
  d = ModiaMedia.density(medium,state)
  h = ModiaMedia.specificEnthalpy(medium,state)

  # Print computed values
  println("data for p=$p, T=$T:")
  println("density          = ", d)
  println("specificEnthalpy = ", h)

  # Plot the most important characteristics of the medium
  ModiaMedia.standardPlot(medium)
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
    constants::AbstractFluidConstants
    limits::FluidLimits
    data  # medium specific data
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
