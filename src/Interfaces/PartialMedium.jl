#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

#=
# A medium is defined by a struct and accompanying functions. The struct has the following structure:
    struct MediumXXX <: AbstractMedium
        mediumName::AbstractString 
        substanceNames::Vector{AbstractString}
        < ... other generic data >
        fluidConstants::Vector{AbstractFluidConstants}
        fluidLimits::FluidLimits
        data  # medium specific data
    end

    All media are stored in dictionary mediumDict
=#



undefinedFunction(name::AbstractString, m::AbstractMedium) = 
    error("\n\nModiaMedia: Function $name not defined for medium $(m.infos.mediumName).\n")

undefinedFunction(name::AbstractString, state::ThermodynamicState) = 
    error("\n\nModiaMedia: Function $name not defined for typeof(state) = ", typeof(state), ".\n")


### Enumerations and constants --------------------------------------------------------------------

"""
    @enum IndependentVariables

Enumeration defining the independent variables of a medium. Possible values:

| Value                      | Independent variables                       |
|:---------------------------|:--------------------------------------------|
| `IndependentVariables_T`   | Temperature                                 |
| `IndependentVariables_pT`  | Pressure, temperature                       |
| `IndependentVariables_ph`  | Pressure, specific enthalpy                 |
| `IndependentVariables_phX` | Pressure, specific enthalpy, mass fractions |
| `IndependentVariables_pTX` | Pressure, temperature, mass fractions       |
| `IndependentVariables_dTX` | Density, temperature, mass fractions        |
"""
@enum IndependentVariables begin 
    IndependentVariables_T=1       # Temperature
    IndependentVariables_pT        # Pressure, Temperature
    IndependentVariables_ph        # Pressure, Specific Enthalpy
    IndependentVariables_phX       # Pressure, Specific Enthalpy, Mass Fraction
    IndependentVariables_pTX       # Pressure, Temperature, Mass Fractions
    IndependentVariables_dTX       # Density, Temperature, Mass Fractions
end



"""
    @enum ReferenceEnthalpy

Enumeration defining the reference enthalpy of a medium. Possible values:

- `ReferenceEnthalpy_ZeroAt0K`: The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded
- `ReferenceEnthalpy_ZeroAt25C`: The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded
- `ReferenceEnthalpy_UserDefined`: The user-defined reference enthalpy is used at 293.15 K (25 degC)
"""
@enum ReferenceEnthalpy begin 
    ReferenceEnthalpy_ZeroAt0K     # The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded
    ReferenceEnthalpy_ZeroAt25C    # The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded
    ReferenceEnthalpy_UserDefined  # The user-defined reference enthalpy is used at 293.15 K (25 degC)
end


"""
    @enum ReferenceEntropy

Enumeration defining the reference entropy of a medium. Possible values:

- `ReferenceEntropy_ZeroAt0K`: The entropy is 0 at 0 K (default)
- `ReferenceEntropy_ZeroAt0C`: The entropy is 0 at 0 degC
- `ReferenceEntropy_UserDefined`: The user-defined reference entropy is used at 293.15 K (25 degC)
"""
@enum ReferenceEntropy begin 
    ReferenceEntropy_ZeroAt0K      # The entropy is 0 at 0 K (default)
    ReferenceEntropy_ZeroAt0C      # The entropy is 0 at 0 degC
    ReferenceEntropy_UserDefined   # The user-defined reference entropy is used at 293.15 K (25 degC)
end


"Enumeration Init defining initialization for fluid flow"
@enum Init begin 
    Init_NoInit                    # NoInit (no initialization)
    Init_InitialStates             # InitialStates (initialize medium states)
    Init_SteadyState               # SteadyState (initialize in steady state)
    Init_SteadyMass                # SteadyMass (initialize density or pressure in steady state)
end


"Enumeration pd defining whether p or d are known for the boundary condition"
@enum pd begin 
    pd_default                     # Default (no boundary condition for p or d)
    pd_p_known                     # p_known (pressure p is known)
    pd_d_known                     # d_known (density d is known)
end


"Enumeration Th - defining whether T or h are known as boundary condition"
@enum Th begin 
    Th_default                     # Default (no boundary condition for T or h)
    Th_T_known                     # T_known (temperature T is known)
    Th_h_known                     # h_known (specific enthalpy h is known)
end




### Data structures ---------------------------------------------------------------------------------

"""
    infos = FluidInfos(; <keyword arguments, see below>)

Generate a new `FluidInfos` object, containing generic properties of a medium.

# Keyword arguments

| Name and type                                  | Description                                                                                             |
|:-----------------------------------------------|:--------------------------------------------------------------------------------------------------------|
| `mediumName::AbstractString`                   | Name of the medium                                                                                      |
| `substanceNames::Vector{AbstractString}`       | Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.                 |
| `extraPropertiesNames::Vector{AbstractString}` | Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused |
| `ThermoStates::IndependentVariables`           | Enumeration type for independent variables                                                              |                                         |
| `singleState::Bool`                            | = true, if u and d are not a function of pressure                                                       |
| `reducedX::Bool`                               | = true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance            |
| `fixedX::Bool`                                 | = true if medium contains the equation X = reference_X                                                  |
| `reference_p`                                  | Reference pressure of Medium: default 1 atmosphere                                                      |
| `reference_T`                                  | Reference temperature of Medium: default 25 deg Celsius                                                 |
| `reference_X`                                  | Reference mass fractions of medium                                                                      |
| `p_default`                                    | Default value for pressure of medium (for initialization)                                               |
| `T_default`                                    | Default value for temperature of medium (for initialization)                                            |
| `h_default`                                    | Default value for specific enthalpy of medium (for initialization)                                      |
| `X_default`                                    | Default value for specific enthalpy of medium (for initialization)                                      |
| `nS`                                           | Number of substances                                                                                    |
| `nX`                                           | Number of mass fractions                                                                                |
| `nXi`                                          | Number of structurally independent mass fractions                                                       |
| `nC`                                           | Number of extra (outside of standard mass-balance) transported properties                               |
| `C_nominal`                                    | Default for the nominal values for the extra properties                                                 |

# Example
```julia
import ModiaMedia

infos = ModiaMedia.FluidInfos(mediumName           = "simpleMedium",
                              substanceNames       = [mediumName],
                              extraPropertiesNames = fill("",0),
                              ThermoStates         = IndependentVariables_T)
```
"""
mutable struct FluidInfos
    mediumName::AbstractString                   # "Name of the medium";
    substanceNames::Vector{AbstractString}       # "Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.";
    extraPropertiesNames::Vector{AbstractString} # "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused"
    ThermoStates::IndependentVariables           # "Enumeration type for independent variables";
    singleState::Bool                            # "= true, if u and d are not a function of pressure";
    reducedX::Bool                               # "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
    fixedX::Bool                                 # "= true if medium contains the equation X = reference_X";
    reference_p::Float64                         # "Reference pressure of Medium: default 1 atmosphere";
    reference_T::Float64                         # "Reference temperature of Medium: default 25 deg Celsius";
    reference_X::AbstractVector                  # "Reference mass fractions of medium";
    p_default::Float64                           # "Default value for pressure of medium (for initialization)";
    T_default::Float64                           # "Default value for temperature of medium (for initialization)";
    h_default::Float64                           # "Default value for specific enthalpy of medium (for initialization)";
    X_default::Vector{Float64}                   # "Default value for specific enthalpy of medium (for initialization)";
    nS::Int                                      # "Number of substances"
    nX::Int                                      # "Number of mass fractions"
    nXi::Int                                     # "Default value for mass fractions of medium (for initialization)"
    nC::Int                                      # "Number of extra (outside of standard mass-balance) transported properties"
    C_nominal::Vector{Float64}                   # "Default for the nominal values for the extra properties"             

    function FluidInfos(; mediumName=Missing,
                          substanceNames=[mediumName],
                          extraPropertiesNames=fill("",0), 
                          ThermoStates=Missing,
                          singleState=Missing, 
                          reducedX=true,
                          fixedX=false, 
                          reference_p=101325.0,
                          reference_T=298.15,
                          reference_X=fill(1.0/length(substanceNames),length(substanceNames)),
                          p_default=101325.0,
                          T_default=293.15,
                          h_default=NaN,
                          X_default=reference_X,
                          C_nominal=1e-6*ones(length(extraPropertiesNames)))
         nS  = length(substanceNames)
         nXi = fixedX ? 0 : ( reducedX ? nS-1 : nS )
         nC  = length(extraPropertiesNames)

         new(mediumName, substanceNames, extraPropertiesNames,
             ThermoStates, singleState, reducedX, fixedX, reference_p,
             reference_T, reference_X, p_default, T_default, h_default, X_default,
             nS, nS, nXi, nC, C_nominal)
    end
end



"""
    fluidLimits = FluidLimits(; TMIN=NaN, TMAX=NaN, DMIN=NaN, DMAX=NaN, PMIN=NaN, PMAX=NaN, 
                                HMIN=NaN, HMAX=NaN, SMIN=NaN, SMAX=NaN)

Generate a new `FluidLimits` object, containing the validity limits of the medium.
"""
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

    FluidLimits(; TMIN=NaN, TMAX=NaN, DMIN=NaN, DMAX=NaN, PMIN=NaN, PMAX=NaN, 
                  HMIN=NaN, HMAX=NaN, SMIN=NaN, SMAX=NaN) =
              new(TMIN, TMAX, DMIN, DMAX, PMIN, PMAX, HMIN, HMAX, SMIN, SMAX)
end



"""
    fluidConstants = BasicFluidConstants(;iupacName="", casRegistryNumber="", 
                         chemicalFormula="", structureFormula="", molarMass=NaN)

Generate a `BasicFluidConstants <: AbstractFluidConstants` object 
containing the minimal information about the standard data of the medium.
"""
mutable struct BasicFluidConstants <: AbstractFluidConstants
    iupacName::AbstractString           # "Complete IUPAC name (or common name, if non-existent)";
    casRegistryNumber::AbstractString   # "Chemical abstracts sequencing number (if it exists)";
    chemicalFormula::AbstractString     # "Chemical formula, (brutto, nomenclature according to Hill";
    structureFormula::AbstractString    # "Chemical structure formula";
    molarMass::Float64                  # "Molar mass";

    BasicFluidConstants(;iupacName="", casRegistryNumber="", chemicalFormula="", structureFormula="", molarMass=NaN) =
       new(iupacName, casRegistryNumber, chemicalFormula, structureFormula, molarMass)
end



"""
    fluidConstants = IdealGasFluidConstants(; iupacName="", casRegistryNumber="", 
          chemicalFormula="", structureFormula="", molarMass=NaN,
          criticalTemperature=NaN, criticalPressure=NaN, criticalMolarVolume=NaN, 
          acentricFactor=NaN, meltingPoint=NaN, normalBoilingPoint=NaN, dipoleMoment=NaN,
          hasIdealGasHeatCapacity=false, hasCriticalData=false, hasDipoleMoment=false,
          hasFundamentalEquation=false, hasLiquidHeatCapacity=false, hasSolidHeatCapacity=false,
          hasAccurateViscosityData=false, hasAccurateConductivityData=false,
          hasVapourPressureCurve=false, hasAcentricFactor=false,
          HCRIT0=0.0, SCRIT0=0.0, deltah=0.0, deltas=0.0)

Generate a `IdealGasFluidConstants <: AbstractFluidConstants` object
containing the minimal information about the standard data of ideal gas media
(critical, triple, molecular and other standard data).
"""
mutable struct IdealGasFluidConstants <: AbstractFluidConstants
    iupacName::AbstractString           # "Complete IUPAC name (or common name, if non-existent)";
    casRegistryNumber::AbstractString   # "Chemical abstracts sequencing number (if it exists)";
    chemicalFormula::AbstractString     # "Chemical formula, (brutto, nomenclature according to Hill";
    structureFormula::AbstractString    # "Chemical structure formula";
    molarMass::Float64                  # "Molar mass";

    criticalTemperature::Float64        # "Critical temperature";
    criticalPressure::Float64           # "Critical pressure";
    criticalMolarVolume::Float64        # "Critical molar Volume";
    acentricFactor::Float64             # "Pitzer acentric factor";
    meltingPoint::Float64               # "Melting point at 101325 Pa";
    normalBoilingPoint::Float64         # "Normal boiling point (at 101325 Pa)";
    dipoleMoment::Float64               # "Dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";

    hasIdealGasHeatCapacity::Bool       # "True if ideal gas heat capacity is available";
    hasCriticalData::Bool               # "True if critical data are known";
    hasDipoleMoment::Bool               # "True if a dipole moment known";
    hasFundamentalEquation::Bool        # "True if a fundamental equation";
    hasLiquidHeatCapacity::Bool         # "True if liquid heat capacity is available";
    hasSolidHeatCapacity::Bool          # "True if solid heat capacity is available";
    hasAccurateViscosityData::Bool      # "True if accurate data for a viscosity function is available";
    hasAccurateConductivityData::Bool   # "True if accurate data for thermal conductivity is available";
    hasVapourPressureCurve::Bool        # "True if vapour pressure data, e.g., Antoine coefficients are known";
    hasAcentricFactor::Bool             # "True if Pitzer accentric factor is known";

    HCRIT0::Float64                     # "Critical specific enthalpy of the fundamental equation";         
    SCRIT0::Float64                     # "Critical specific entropy of the fundamental equation";
    deltah::Float64                     # "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
    deltas::Float64                     # "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";

    IdealGasFluidConstants(; iupacName="", casRegistryNumber="", chemicalFormula="", structureFormula="", molarMass=NaN,
                             criticalTemperature=NaN, criticalPressure=NaN, criticalMolarVolume=NaN, 
                             acentricFactor=NaN, meltingPoint=NaN, normalBoilingPoint=NaN, dipoleMoment=NaN,
                             hasIdealGasHeatCapacity=false, hasCriticalData=false, hasDipoleMoment=false,
                             hasFundamentalEquation=false, hasLiquidHeatCapacity=false, hasSolidHeatCapacity=false,
                             hasAccurateViscosityData=false, hasAccurateConductivityData=false,
                             hasVapourPressureCurve=false, hasAcentricFactor=false,
                             HCRIT0=0.0, SCRIT0=0.0, deltah=0.0, deltas=0.0) =
        new(iupacName, casRegistryNumber, chemicalFormula, structureFormula, molarMass,
            criticalTemperature, criticalPressure, criticalMolarVolume, acentricFactor, meltingPoint,
            normalBoilingPoint, dipoleMoment, hasIdealGasHeatCapacity, hasCriticalData, hasDipoleMoment,
            hasFundamentalEquation, hasLiquidHeatCapacity, hasSolidHeatCapacity, hasAccurateViscosityData,
            hasAccurateConductivityData, hasVapourPressureCurve, hasAcentricFactor,
            HCRIT0, SCRIT0, deltah, deltas)
end


""" 
    state = setState_pTX(medium, p,T,X)

Generate a state object for medium `medium::AbstractMedium` for
pressure `p` [Pa], temperature `T` [K] and mass fractions vector `X` or `Xi`.
"""
setState_pTX(m::AbstractMedium,p,T,X) = undefinedFunction("setState_pTX", m)



""" 
    setState_pTX!(state, p,T,X)

Update the `state::ThermodynamicState` object with 
pressure `p` [Pa], temperature `T` [K] and mass fractions vector `X` or `Xi`.
"""
setState_pTX!(state::ThermodynamicState,p,T,X) = undefinedFunction("setState_pTX!", state)



""" 
    state = setState_phX(medium, p,h,X)

Generate a state object for medium `medium::AbstractMedium` for
pressure `p` [Pa], specific enthalpy `h` [J/kg]] and mass fractions vector `X` or `Xi`.
"""
setState_phX(m::AbstractMedium,p,h,X) = undefinedFunction("setState_phX", m)



""" 
    setState_phX!(state, p,h,X)

Update the `state::ThermodynamicState` object with 
pressure `p` [Pa], specific enthalpy `h` [J/kg]] and mass fractions vector `X` or `Xi`.
"""
setState_phX!(state::ThermodynamicState,p,h,X) = undefinedFunction("setState_phX!", state)




""" 
    state = setState_psX(medium, p,s,X)

Generate a state object for medium `medium::AbstractMedium` for
pressure `p` [Pa], specific entropy `s` [J/(kg*K)] and mass fractions vector `X`  or `Xi`.
"""
setState_psX(m::AbstractMedium,p,s,X) = undefinedFunction("setState_psX", m)



""" 
    setState_psX!(state, p,s,X)

Update the `state::ThermodynamicState` object with 
pressure `p` [Pa], specific entropy `s` [J/(kg*K)] and mass fractions vector `X` or `Xi`.
"""
setState_psX!(state::ThermodynamicState,p,s,X) = undefinedFunction("setState_psX!", state)



""" 
    state = setState_dTX(medium, d,T,X)

Generate a state object for medium `medium::AbstractMedium` for
density `d` [kg/m^3], temperature `T` [K] and mass fractions vector `X` or `Xi`.
"""
setState_dTX(m::AbstractMedium,d,T,X) = undefinedFunction("setState_dTX", m)



""" 
    setState_dTX!(state, d,T,X)

Update the `state::ThermodynamicState` object with 
density `d` [kg/m^3], temperature `T` [K] and mass fractions vector `X` or `Xi`.
"""
setState_dTX!(state::ThermodynamicState,d,T,X) = undefinedFunction("setState_dTX!", state)



""" 
    state_b = isenthalpicState(state_a,dp)

Return `state_b` by an isenthalpic transformation of `state_a` with pressure drop dp:

```julia
pressure(state_b)         = pressure(state_a) + dp
specificEnthalpy(state_b) = specificEnthalpy(state_a)
state_b.X                 = state_a.X
```
"""
isenthalpicState(m::AbstractMedium, state_a::ThermodynamicState, dp::Float64) = undefinedFunction("isenthalpicState", m)
isenthalpicState(                   state_a::ThermodynamicState, dp::Float64) = isenthalpicState(state_a.Medium, state_a, dp)


""" 
    isenthalpicState!(state_b,state_a,dp)

Update `state_b` by an isenthalpic transformation of `state_a` with pressure drop dp:

```julia
pressure(state_b)         = pressure(state_a) + dp
specificEnthalpy(state_b) = specificEnthalpy(state_a)
state_b.X                 = state_a.X
```
"""
isenthalpicState!(state_b::ThermodynamicState, state_a::ThermodynamicState, dp::Float64) = undefinedFunction("isenthalpicState!", state)


"""
    getMedium(state)

Return the Medium of thermodynamic `state`
"""
getMedium(state::ThermodynamicState) = state.Medium



### Medium functions ---------------------------------------------------------------------------------

"""
    pressure(state)

Return pressure from `state::ThermodynamicState` in [Pa]
"""
pressure(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("pressure", m)
pressure(                   state::ThermodynamicState)::Float64 = pressure(state.Medium, state)


"""
    temperature(state)

Return temperature from `state::ThermodynamicState` in [K]
"""
temperature(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("temperature", m)
temperature(                   state::ThermodynamicState)::Float64 = temperature(state.Medium, state)



""" 
    density(state)

Return density from `state::ThermodynamicState` in [kg/m^3]
"""
density(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("density", m)
density(                   state::ThermodynamicState)::Float64 = density(state.Medium, state)



"""
    specificEnthalpy(state)

Return specific enthalpy from `state::ThermodynamicState` in [J/kg]
"""
specificEnthalpy(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("specificEnthalpy", m)
specificEnthalpy(                   state::ThermodynamicState)::Float64 = specificEnthalpy(state.Medium, state)



"""
    specificInternalEnergy(state)

Return specific internal energy from `state::ThermodynamicState` in [J/kg]
"""
specificInternalEnergy(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("specificInternalEnergy", m)
specificInternalEnergy(                   state::ThermodynamicState)::Float64 = specificInternalEnergy(state.Medium, state)



""" 
    specificHeatCapacityCp(state)

Return specific heat capacity at constant pressure from `state::ThermodynamicState` in [J/(kg*K)]
"""
specificHeatCapacityCp(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("specificHeatCapacity", m)
specificHeatCapacityCp(                   state::ThermodynamicState)::Float64 = specificHeatCapacityCp(state.Medium, state)


"""
    dynamicViscosity(state)

Return dynamic viscosity from `state::ThermodynamicState` in [Pa*s]
"""
dynamicViscosity(m::AbstractMedium, state::ThermodynamicState)::Float64 = undefinedFunction("dynamicViscosity", m)
dynamicViscosity(                   state::ThermodynamicState)::Float64 = dynamicViscosity(state.Medium, state)





### Medium functions expressed from the functions above -------------------------------------------------

"""
    specificEnthalpy_pTX(medium,p,T,X)

Return specific enthalpy for `medium::AbstractMedium` from p, T, and X or Xi in [J/kg].
"""
specificEnthalpy_pTX(m::AbstractMedium, p,T,X) = specficEnthalpy(setState_pTX(m,p,T,X))



""" 
    temperature_phX(medium,p,h,X)

Return temperature for `medium::AbstractMedium` from p, h, and X or Xi in [K].
"""
temperature_phX(m::AbstractMedium, p,h,X) = temperature(setState_phX(m,p,h,X))


"""
    density_pTX(medium,p,T,X)

Return density for `medium::AbstractMedium` from p, T, and X or Xi in [kg/m^3]
"""
density_pTX(m::AbstractMedium, p,T,X) = density(setState_pTX(m,p,T,X))



"""
    density_phX(medium,p,h,X)

Return density for `medium::AbstractMedium` from p, h, and X or Xi in [kg/m^3]
"""
density_phX(m::AbstractMedium, p,h,X) = density(setState_phX(m,p,h,X))



### Other functions operating on the medium ---------------------------------------------------------------------------------

to_DensityDisplayUnit(d) = d*1e-3
from_degC(Celsius)       = Celsius + 273.15


"""
    dict = standardCharacteristics(medium::AbstractMedium)

Return a `dict::dict{AbstractString,Any}` dictionary with the most 
important characteristics of the medium as vectors.
"""
standardCharacteristics(m::AbstractMedium)::dict{AbstractString,Any} = undefinedFunction("standardCharacteristics", m)


"""
    standardPlot(medium::AbstractMedium; figure=1)

Plot the `standardCharacteristics(medium)` of the medium.
"""
standardPlot(m::AbstractMedium; figure=1) = undefinedFunction("standardPlot", m)


const T_zero = -273.15
to_degC(T_in_K) = T_in_K + Modelica.Constants.T_zero;




