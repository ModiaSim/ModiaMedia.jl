# ModiaMedia

[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://modiasim.github.io/ModiaMedia.jl/latest/)

This package  provides thermodynamic property models for use with [Modia](https://github.com/ModiaSim/Modia.jl)
and other Julia packages. The initial goal is to achieve a similar functionality as
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media),
the standard media library for Modelica models, but with improvements based on Julia features
such as multiple dispatch.

This package is under development and it is planned to provide all thermodynamic property models from
[Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media)
in this package.


## Installation

This package is currently under development and is not yet registered in METADATA.
Julia 1.7 is required. Installation is performed via:

```julia
julia> ]add https://github.com/ModiaSim/ModiaMedia.jl
```

ModiaMedia performs plotting via [SignalTables.jl](https://github.com/ModiaSim/SignalTables.jl).

## Use

```julia
  using ModiaMedia

  # Define thermodynamic property model to be used
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
  using SignalTables
  usePlotPackage("PyPlot")   # for details see SignalTables
  @usingPlotPackage          # using SignalTablesInterface_PyPlot
  standardPlot(Medium, plot)
```

The last command results in the following plot:

![standardPlot](https://ModiaSim.github.io/ModiaMedia.jl/resources/images/N2.png)


## Main Developers

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/) ([DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en))
- Hilding Elmqvist ([Mogram](http://www.mogram.net/)),
- [Chris Laughman](http://www.merl.com/people/laughman) ([MERL](http://www.merl.com/)).
- All the content of ModiaMedia is based on
  [Modelica.Media](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Media.html#Modelica.Media)
  which was and is developed from many people.

License: MIT (expat)
