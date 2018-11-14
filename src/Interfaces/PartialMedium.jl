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
    IndependentVariables_T         # Temperature
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
| `ThermoStates::IndependentVariables`           | Enumeration type for independent variables                                                              |
| `baseProperties::Symbol`                       | Symbol of baseProperties model = :BaseProperties_<StructName>                                           |
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
                          baseProperties=Missing,
                          singleState=Missing, 
                          reducedX=true,
                          fixedX=false, 
                          reference_p=101325,
                          reference_T=298.15,
                          reference_X=fill(1/length(substanceNames),length(substanceNames)),
                          p_default=101325,
                          T_default=293.15,
                          h_default=NaN,
                          X_default=reference_X,
                          C_nominal=1e-6*ones(length(extraPropertiesNames)))
         nS  = length(substanceNames)
         nXi = fixedX ? 0 : ( reducedX ? nS-1 : nS )
         nC  = length(extraPropertiesNames)

         new(mediumName, substanceNames, extraPropertiesNames,
             ThermoStates, baseProperties, singleState, reducedX, fixedX, reference_p,
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
struct BasicFluidConstants <: AbstractFluidConstants
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
    state = setState_phX(medium, p,h,X)

Generate a state object for medium `medium::AbstractMedium` for
pressure `p` [Pa], specific enthalpy `h` [J/kg]] and mass fractions vector `X` or `Xi`.
"""
setState_phX(m::AbstractMedium,p,h,X) = undefinedFunction("setState_phX", m)



""" 
    state = setState_psX(medium, p,s,X)

Generate a state object for medium `medium::AbstractMedium` for
pressure `p` [Pa], specific entropy `s` [J/(kg*K)] and mass fractions vector `X`  or `Xi`.
"""
setState_psX(m::AbstractMedium,p,s,X) = undefinedFunction("setState_psX", m)


""" 
    state = setState_dTX(medium, d,T,X)

Generate a state object for medium `medium::AbstractMedium` for
density `d` [kg/m^3], temperature `T` [K] and mass fractions vector `X` or `Xi`.
"""
setState_dTX(m::AbstractMedium,d,T,X) = undefinedFunction("setState_dTX", m)



### Medium functions ---------------------------------------------------------------------------------

"p = pressure(medium,state) - return pressure for `medium` from `state` in [Pa]"
pressure(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("pressure", m)

"T = temperature(medium,state) - return temperature for `medium` from `state` in [K]"
temperature(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("temperature", m)

"d = density(medium,state) - return density for `medium` from `state` in [kg/m^3]"
density(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("density", m)

"h = specificEnthalpy(medium,state) - return specific enthalpy for `medium` from `state` in [J/kg]"
specificEnthalpy(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("specificEnthalpy", m)

"u = specificInternalEnergy(medium,state) - return specific internal energy from `medium` at `state` in [J/kg]"
specificInternalEnergy(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("specificInternalEnergy", m)

"cp = specificHeatCapacityCp(medium,state) - return specific heat capacity at constant pressure for `medium` from `state` in [J/(kg*K)]"
specificHeatCapacityCp(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("specificHeatCapacity", m)



### Medium functions expressed from the functions above -------------------------------------------------

"h = specificEnthalpy_pTX(medium,p,T,X) - return specific enthalpy for `medium` from p, T, and X or Xi in [J/kg]"
specificEnthalpy_pTX(m::AbstractMedium, p,T,X) = specficEnthalpy(m, setState_pTX(m,p,T,X))

"T = temperature_phX(medium,p,h,X) - return temperature for `medium` from p, h, and X or Xi in [K]"
temperature_phX(m::AbstractMedium, p,h,X) = temperature(m, setState_phX(m,p,h,X))

"d = density_pTX(medium,p,T,X) - return density for `medium` from p, T, and X or Xi in [kg/m^3]"
density_pTX(m::AbstractMedium, p,T,X) = density(m, setState_pTX(m,p,T,X))

"d = density_phX(medium,p,h,X) - return density for `medium` from p, h, and X or Xi in [kg/m^3]"
density_phX(m::AbstractMedium, p,h,X) = density(m, setState_phX(m,p,h,X))



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




