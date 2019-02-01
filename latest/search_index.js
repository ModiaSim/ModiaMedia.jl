var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ModiaMedia.jl-Documentation-1",
    "page": "Home",
    "title": "ModiaMedia.jl Documentation",
    "category": "section",
    "text": "ModiaMedia shall provide Media models for use with Modia and other Julia packages. The initial goal is to achieve a similar functionality as Modelica.Media, the standard media library for Modelica models, but with improvements based on Julia features such as multiple dispatch.This package is under development and it is planned to provide all media from Modelica.Media in this package."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package is currently under development and is not yet registered in METADATA. Julia 1.0 is required. Installation is performed via:julia> ]add https://github.com/ModiaSim/ModiaMedia.jlModiaMedia uses PyPlot for plotting (via ModiaMath.plot). If PyPlot is not available in your current Julia environment an information message is printed and all plot(..) calls are ignored.In order that plot windows are displayed, you need to add PyPlot to your current environment via ]add PyPlot. Often this automatic installation fails and it is recommended to follow the instructions Installing PyPlot in a robust way."
},

{
    "location": "index.html#Use-1",
    "page": "Home",
    "title": "Use",
    "category": "section",
    "text": "  using ModiaMedia\r\n\r\n  # Define medium to be used\r\n  Medium = getMedium(\"N2\");\r\n\r\n  # Define the operating point where the medium shall be evaluated.\r\n  p = 1e5    # in [Pa]\r\n  T = 300.0  # in [K]\r\n\r\n  # Set the medium-specific thermodynamic state from p and T\r\n  # (could be also set from p and h, or p and s, or d and T)\r\n  state = setState_pT(Medium, p, T)\r\n\r\n  # Call media functions (here to compute density and specific enthalpy)\r\n  d = density(Medium,state)\r\n  h = specificEnthalpy(Medium,state)\r\n\r\n  # Print computed values\r\n  println(\"data for p=$p, T=$T:\")\r\n  println(\"density          = \", d)\r\n  println(\"specificEnthalpy = \", h)\r\n\r\n  # Plot the most important characteristics of the medium\r\n  ModiaMedia.standardPlot(Medium)This example generates the following plot:(Image: standardPlot)"
},

{
    "location": "index.html#Available-media-1",
    "page": "Home",
    "title": "Available media",
    "category": "section",
    "text": "The available media can be listed with listMedia() resulting in:Media available in ModiaMedia:\r\n\r\n¦ Row ¦ name                        ¦ type                 ¦\r\n+-----+-----------------------------+----------------------¦\r\n¦ 1   ¦ Ar                          ¦ SingleGasNasa        ¦\r\n¦ 2   ¦ C2H2_vinylidene             ¦ SingleGasNasa        ¦\r\n¦ 3   ¦ C2H4                        ¦ SingleGasNasa        ¦\r\n¦ 4   ¦ C2H5OH                      ¦ SingleGasNasa        ¦\r\n¦ 5   ¦ C2H6                        ¦ SingleGasNasa        ¦\r\n¦ 6   ¦ C3H6_propylene              ¦ SingleGasNasa        ¦\r\n¦ 7   ¦ C3H8                        ¦ SingleGasNasa        ¦\r\n¦ 8   ¦ C4H10_n_butane              ¦ SingleGasNasa        ¦\r\n¦ 9   ¦ C4H8_1_butene               ¦ SingleGasNasa        ¦\r\n¦ 10  ¦ C5H10_1_pentene             ¦ SingleGasNasa        ¦\r\n¦ 11  ¦ C5H12_n_pentane             ¦ SingleGasNasa        ¦\r\n¦ 12  ¦ C6H12_1_hexene              ¦ SingleGasNasa        ¦\r\n¦ 13  ¦ C6H14_n_hexane              ¦ SingleGasNasa        ¦\r\n¦ 14  ¦ C6H6                        ¦ SingleGasNasa        ¦\r\n¦ 15  ¦ C7H14_1_heptene             ¦ SingleGasNasa        ¦\r\n¦ 16  ¦ C7H16_n_heptane             ¦ SingleGasNasa        ¦\r\n¦ 17  ¦ C8H10_ethylbenz             ¦ SingleGasNasa        ¦\r\n¦ 18  ¦ C8H18_n_octane              ¦ SingleGasNasa        ¦\r\n¦ 19  ¦ CH3OH                       ¦ SingleGasNasa        ¦\r\n¦ 20  ¦ CH4                         ¦ SingleGasNasa        ¦\r\n¦ 21  ¦ CL2                         ¦ SingleGasNasa        ¦\r\n¦ 22  ¦ CO                          ¦ SingleGasNasa        ¦\r\n¦ 23  ¦ CO2                         ¦ SingleGasNasa        ¦\r\n¦ 24  ¦ ConstantPropertyLiquidWater ¦ SimpleMedium         ¦\r\n¦ 25  ¦ F2                          ¦ SingleGasNasa        ¦\r\n¦ 26  ¦ H2                          ¦ SingleGasNasa        ¦\r\n¦ 27  ¦ H2O                         ¦ SingleGasNasa        ¦\r\n¦ 28  ¦ He                          ¦ SingleGasNasa        ¦\r\n¦ 29  ¦ MoistAir                    ¦ MoistAir             ¦\r\n¦ 30  ¦ N2                          ¦ SingleGasNasa        ¦\r\n¦ 31  ¦ N2O                         ¦ SingleGasNasa        ¦\r\n¦ 32  ¦ NH3                         ¦ SingleGasNasa        ¦\r\n¦ 33  ¦ NO                          ¦ SingleGasNasa        ¦\r\n¦ 34  ¦ NO2                         ¦ SingleGasNasa        ¦\r\n¦ 35  ¦ Ne                          ¦ SingleGasNasa        ¦\r\n¦ 36  ¦ O2                          ¦ SingleGasNasa        ¦\r\n¦ 37  ¦ SO2                         ¦ SingleGasNasa        ¦\r\n¦ 38  ¦ SO3                         ¦ SingleGasNasa        ¦\r\n¦ 39  ¦ SimpleAir                   ¦ SimpleIdealGasMedium ¦"
},

{
    "location": "index.html#Structure-of-package-1",
    "page": "Home",
    "title": "Structure of package",
    "category": "section",
    "text": "A medium is a struct of the following type:mutable struct MediumName <: ModiaMedia.AbstractMedium  # or of a subtype of AbstractMedium\r\n    infos::ModiaMedia.FluidInfos\r\n    fluidConstants::Vector{ModiaMedia.AbstractFluidConstants}\r\n    fluidLimits::ModiaMedia.FluidLimits\r\n    data  # medium specific data\r\nend\r\n\r\nstruct FluidInfos\r\n    mediumName::AbstractString                   # \"Name of the medium\";\r\n    substanceNames::Vector{AbstractString}       # \"Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.\";\r\n    extraPropertiesNames::Vector{AbstractString} # \"Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\\\"\\\",0) if unused\"\r\n    ThermoStates::IndependentVariables           # \"Enumeration type for independent variables\";\r\n    singleState::Bool                            # \"= true, if u and d are not a function of pressure\";\r\n    reducedX::Bool                               # \"= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)\";\r\n    fixedX::Bool                                 # \"= true if medium contains the equation X = reference_X\";\r\n    reference_p::Float64                         # \"Reference pressure of Medium: default 1 atmosphere\";\r\n    reference_T::Float64                         # \"Reference temperature of Medium: default 25 deg Celsius\";\r\n    reference_X::AbstractVector                  # \"Default mass fractions of medium\";\r\n    p_default::Float64                           # \"Default value for pressure of medium (for initialization)\";\r\n    T_default::Float64                           # \"Default value for temperature of medium (for initialization)\";\r\n    h_default::Float64                           # \"Default value for specific enthalpy of medium (for initialization)\";\r\n    X_default::Vector{Float64}                   # \"Default value for specific enthalpy of medium (for initialization)\";\r\n    nS::Int                                      # \"Number of substances\"\r\n    nX::Int                                      # \"Number of mass fractions\"\r\n    nXi::Int                                     # \"Default value for mass fractions of medium (for initialization)\"\r\n    nC::Int                                      # \"Number of extra (outside of standard mass-balance) transported properties\"\r\n    C_nominal::Vector{Float64}                   # \"Default for the nominal values for the extra properties\"\r\nend\r\n\r\nstruct BasicFluidConstants <: AbstractFluidConstants\r\n    iupacName::AbstractString           # \"Complete IUPAC name (or common name, if non-existent)\";\r\n    casRegistryNumber::AbstractString   # \"Chemical abstracts sequencing number (if it exists)\";\r\n    chemicalFormula::AbstractString     # \"Chemical formula, (brutto, nomenclature according to Hill\";\r\n    structureFormula::AbstractString    # \"Chemical structure formula\";\r\n    molarMass::Float64                  # \"Molar mass\";\r\nend\r\n\r\nmutable struct FluidLimits\r\n    TMIN::Float64  # \"Minimum temperature\";\r\n    TMAX::Float64  # \"Maximum temperature\";\r\n    DMIN::Float64  # \"Minimum density\";\r\n    DMAX::Float64  # \"Maximum density\";\r\n    PMIN::Float64  # \"Minimum pressure\";\r\n    PMAX::Float64  # \"Maximum pressure\";\r\n    HMIN::Float64  # \"Minimum enthalpy\";\r\n    HMAX::Float64  # \"Maximum enthalpy\";\r\n    SMIN::Float64  # \"Minimum entropy\";\r\n    SMAX::Float64  # \"Maximum entropy\";\r\nendand all instances of this struct are stored in a dictionary. This dictionary is constructed in a preprocessing step by running \"ModiaMedia/dict/GenerateMediumDict.jl\". This module contains code that was mostly automatically converted from Modelica.Media to Julia. The resulting dictionary is serialized and stored in \"ModiaMedia/src/Media/media.julia_serializer\". When package ModiaMedia is compiled, this serialized dictionary is deserialized and included in the compiled package.Function ModiaMedia.Medium(name) returns the MediumXXX instance stored in the medium dictionary with key name."
},

{
    "location": "index.html#Main-Developers-1",
    "page": "Home",
    "title": "Main Developers",
    "category": "section",
    "text": "Martin Otter (DLR - Institute of System Dynamics and Control)\\\nHilding Elmqvist (Mogram),\\\nChris Laughman (MERL).License: MIT (expat)"
},

{
    "location": "index.html#Release-Notes-1",
    "page": "Home",
    "title": "Release Notes",
    "category": "section",
    "text": "The ModiaMedia package development has just started and a lot has to be improved."
},

{
    "location": "index.html#Version-0.1.0-dev-1",
    "page": "Home",
    "title": "Version 0.1.0-dev",
    "category": "section",
    "text": "A version is not yet released."
},

{
    "location": "lib/Types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Types.html#ModiaMedia.AbstractMedium",
    "page": "Types",
    "title": "ModiaMedia.AbstractMedium",
    "category": "type",
    "text": "abstract type AbstractMedium - Abstract type of all media\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IndependentVariables",
    "page": "Types",
    "title": "ModiaMedia.IndependentVariables",
    "category": "type",
    "text": "@enum IndependentVariables\n\nEnumeration defining the independent variables of a medium. Possible values:\n\nValue Independent variables\nIndependentVariables_T Temperature\nIndependentVariables_pT Pressure, temperature\nIndependentVariables_ph Pressure, specific enthalpy\nIndependentVariables_phX Pressure, specific enthalpy, mass fractions\nIndependentVariables_pTX Pressure, temperature, mass fractions\nIndependentVariables_dTX Density, temperature, mass fractions\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MoistAir",
    "page": "Types",
    "title": "ModiaMedia.MoistAir",
    "category": "type",
    "text": "medium = MoistAir()\n\nGenerate an MoistAir <: CondensingGases medium object. Valid for T = 190 K ... 647 K. \n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.PureSubstance",
    "page": "Types",
    "title": "ModiaMedia.PureSubstance",
    "category": "type",
    "text": "abstract type PureSubstance <: AbstractMedium - Abstract type of all media consisting of a pure substance\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleIdealGasMedium",
    "page": "Types",
    "title": "ModiaMedia.SimpleIdealGasMedium",
    "category": "type",
    "text": "medium = SimpleIdealGasMedium(; mediumName     = nothing,\n                                reference_p    = 101325.0,\n                                reference_T    = 298.15,\n                                p_default      = 101325.0,\n                                T_default      = 293.15,\n                                fluidConstants = nothing, \n                                data           = nothing)\n\nGenerate a SimpleIdealGasMedium <: PureSubstance medium object.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleMedium",
    "page": "Types",
    "title": "ModiaMedia.SimpleMedium",
    "category": "type",
    "text": "medium = SimpleMedium(; mediumName     = nothing,\n                        reference_p    = 101325.0,\n                        reference_T    = 298.15,\n                        p_default      = 101325.0,\n                        T_default      = 293.15,\n                        fluidConstants = nothing, \n                        data           = nothing)\n\nGenerate a SimpleMedium <: PureSubstance medium object.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasa",
    "page": "Types",
    "title": "ModiaMedia.SingleGasNasa",
    "category": "type",
    "text": "medium = SingleGasNasa(; mediumName     = Missing,\n                         reference_p    = 101325,\n                         reference_T    = 298.15,\n                         p_default      = 101325,\n                         T_default      = 293.15,\n                         fluidConstants = nothing, \n                         fluidLimits    = FluidLimits(TMIN=200.0, TMAX=6000.0), \n                         data           = nothing)\n\nGenerate a SingleGasNasa <: PureSubstance medium object.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ThermodynamicState",
    "page": "Types",
    "title": "ModiaMedia.ThermodynamicState",
    "category": "type",
    "text": "abstract type ThermodynamicState - Abstract type of all media states\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ThermodynamicState_pT",
    "page": "Types",
    "title": "ModiaMedia.ThermodynamicState_pT",
    "category": "type",
    "text": "state = ThermodynamicState_pT(Medium,p,T)\n\nGenerate a ThermodynamicState_pT <: ThermodynamicState object containg pressure p [Pa] and temperature T [K] as states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.AbstractFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.AbstractFluidConstants",
    "category": "type",
    "text": "abstract type AbstractFluidConstants - Abstract type of all FluidConstants structures\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.BasicFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.BasicFluidConstants",
    "category": "type",
    "text": "fluidConstants = BasicFluidConstants(;iupacName=\"\", casRegistryNumber=\"\", \n                     chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN)\n\nGenerate a BasicFluidConstants <: AbstractFluidConstants object  containing the minimal information about the standard data of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.CondensingGases",
    "page": "Types",
    "title": "ModiaMedia.CondensingGases",
    "category": "type",
    "text": "abstract type CondensingGases <: AbstractMedium - Abstract type of all media consisting of condensing media\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidInfos",
    "page": "Types",
    "title": "ModiaMedia.FluidInfos",
    "category": "type",
    "text": "infos = FluidInfos(; <keyword arguments, see below>)\n\nGenerate a new FluidInfos object, containing generic properties of a medium.\n\nKeyword arguments\n\nName and type Description\nmediumName::AbstractString Name of the medium\nsubstanceNames::Vector{AbstractString} Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.\nextraPropertiesNames::Vector{AbstractString} Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused\nThermoStates::IndependentVariables Enumeration type for independent variables\nsingleState::Bool = true, if u and d are not a function of pressure\nreducedX::Bool = true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance\nfixedX::Bool = true if medium contains the equation X = reference_X\nreference_p Reference pressure of Medium: default 1 atmosphere\nreference_T Reference temperature of Medium: default 25 deg Celsius\nreference_X Reference mass fractions of medium\np_default Default value for pressure of medium (for initialization)\nT_default Default value for temperature of medium (for initialization)\nh_default Default value for specific enthalpy of medium (for initialization)\nX_default Default value for specific enthalpy of medium (for initialization)\nnS Number of substances\nnX Number of mass fractions\nnXi Number of structurally independent mass fractions\nnC Number of extra (outside of standard mass-balance) transported properties\nC_nominal Default for the nominal values for the extra properties\n\nExample\n\nimport ModiaMedia\n\ninfos = ModiaMedia.FluidInfos(mediumName           = \"simpleMedium\",\n                              substanceNames       = [mediumName],\n                              extraPropertiesNames = fill(\"\",0),\n                              ThermoStates         = IndependentVariables_T)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidLimits",
    "page": "Types",
    "title": "ModiaMedia.FluidLimits",
    "category": "type",
    "text": "fluidLimits = FluidLimits(; TMIN=NaN, TMAX=NaN, DMIN=NaN, DMAX=NaN, PMIN=NaN, PMAX=NaN, \n                            HMIN=NaN, HMAX=NaN, SMIN=NaN, SMAX=NaN)\n\nGenerate a new FluidLimits object, containing the validity limits of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IdealGasFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.IdealGasFluidConstants",
    "category": "type",
    "text": "fluidConstants = IdealGasFluidConstants(; iupacName=\"\", casRegistryNumber=\"\", \n      chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN,\n      criticalTemperature=NaN, criticalPressure=NaN, criticalMolarVolume=NaN, \n      acentricFactor=NaN, meltingPoint=NaN, normalBoilingPoint=NaN, dipoleMoment=NaN,\n      hasIdealGasHeatCapacity=false, hasCriticalData=false, hasDipoleMoment=false,\n      hasFundamentalEquation=false, hasLiquidHeatCapacity=false, hasSolidHeatCapacity=false,\n      hasAccurateViscosityData=false, hasAccurateConductivityData=false,\n      hasVapourPressureCurve=false, hasAcentricFactor=false,\n      HCRIT0=0.0, SCRIT0=0.0, deltah=0.0, deltas=0.0)\n\nGenerate a IdealGasFluidConstants <: AbstractFluidConstants object containing the minimal information about the standard data of ideal gas media (critical, triple, molecular and other standard data).\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.Init",
    "page": "Types",
    "title": "ModiaMedia.Init",
    "category": "type",
    "text": "Enumeration Init defining initialization for fluid flow\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MixtureMedium",
    "page": "Types",
    "title": "ModiaMedia.MixtureMedium",
    "category": "type",
    "text": "abstract type MixtureMedium <: AbstractMedium - Abstract type of all media consisting of a mixture\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MixtureThermodynamicState",
    "page": "Types",
    "title": "ModiaMedia.MixtureThermodynamicState",
    "category": "type",
    "text": "abstract type MixtureThermodynamicState <: ThermodynamicState - Abstract type of the states of all media consisting of a mixture\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MoistAirData",
    "page": "Types",
    "title": "ModiaMedia.MoistAirData",
    "category": "type",
    "text": "data = MoistAirData(steam::SingleGasNasaData, dryair::SingleGasNasaData)\n\nGenerate an MoistAirData object containing the data for a moist air model.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MoistAirState",
    "page": "Types",
    "title": "ModiaMedia.MoistAirState",
    "category": "type",
    "text": "state = MoistAirState(Medium, p, T, X)\n\nGenerate an MoistAirState <: MixtureThermodynamicState object containing pressure p [Pa], temperature T [K], and a vector of mass fractions as states, where X[1] is the mass fraction of Steam and X[2] is the mass fraction of dry air. If argument X has only one element, X[2] is computed from X[1] and stored in the state.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ReferenceEnthalpy",
    "page": "Types",
    "title": "ModiaMedia.ReferenceEnthalpy",
    "category": "type",
    "text": "@enum ReferenceEnthalpy\n\nEnumeration defining the reference enthalpy of a medium. Possible values:\n\nReferenceEnthalpy_ZeroAt0K: The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded\nReferenceEnthalpy_ZeroAt25C: The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded\nReferenceEnthalpy_UserDefined: The user-defined reference enthalpy is used at 293.15 K (25 degC)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ReferenceEntropy",
    "page": "Types",
    "title": "ModiaMedia.ReferenceEntropy",
    "category": "type",
    "text": "@enum ReferenceEntropy\n\nEnumeration defining the reference entropy of a medium. Possible values:\n\nReferenceEntropy_ZeroAt0K: The entropy is 0 at 0 K (default)\nReferenceEntropy_ZeroAt0C: The entropy is 0 at 0 degC\nReferenceEntropy_UserDefined: The user-defined reference entropy is used at 293.15 K (25 degC)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleIdealGasMediumData",
    "page": "Types",
    "title": "ModiaMedia.SimpleIdealGasMediumData",
    "category": "type",
    "text": "data = SimpleIdealGasMediumData(;cp_const=nothing, R_gas=nothing, MM_const=nothing, \n                                 eta_const=nothing, lambda_const=nothing, T_min=nothing,\n                                 T_max=nothing, T0=nothing)\n\nGenerate a SimpleIdealGasMediumData object containing the data for a SimpleIdealGas medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleIdealGasMediumState",
    "page": "Types",
    "title": "ModiaMedia.SimpleIdealGasMediumState",
    "category": "type",
    "text": "state = SimpleIdealGasMediumState(Medium, p, T)\n\nGenerate a SimpleIdealGasMediumState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleMediumData",
    "page": "Types",
    "title": "ModiaMedia.SimpleMediumData",
    "category": "type",
    "text": "data = SimpleMediumData(;cp_const=nothing, cv_const=nothing, d_const=nothing, eta_const=nothing,\n                         lambda_const=nothing, a_const=nothing, T0=nothing, MM_const=nothing)\n\nGenerate a SimpleMediumData object containing the data for a SimpleMedium medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleMediumState",
    "page": "Types",
    "title": "ModiaMedia.SimpleMediumState",
    "category": "type",
    "text": "state = SimpleMediumState(Medium, p, T)\n\nGenerate a SimpleMediumState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasaData",
    "page": "Types",
    "title": "ModiaMedia.SingleGasNasaData",
    "category": "type",
    "text": "data = SingleGasNasaData(;name=Missing, MM=nothing, Hf=nothing, H0=nothing, Tlimit=nothing,\n                          alow=Missing, blow=Missing, ahigh=Missing,\n                          bhigh=Missing, R=nothing)\n\nGenerate a SingleGasNasaData object containing the data for an ideal Gas based on the NASA Glenn coefficients.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasaState",
    "page": "Types",
    "title": "ModiaMedia.SingleGasNasaState",
    "category": "type",
    "text": "state = SingleGasNasaState(Medium, p, T)\n\nGenerate a SingleGasNasaState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.Th",
    "page": "Types",
    "title": "ModiaMedia.Th",
    "category": "type",
    "text": "Enumeration Th - defining whether T or h are known as boundary condition\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ThermodynamicState_pTX",
    "page": "Types",
    "title": "ModiaMedia.ThermodynamicState_pTX",
    "category": "type",
    "text": "state = ThermodynamicState_pTX(p,T,X)\n\nGenerate a ThermodynamicState_pT <: MixtureThermodynamicState object containing pressure p [Pa], temperature T [K], and a vector of mass fractions X as states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.pd",
    "page": "Types",
    "title": "ModiaMedia.pd",
    "category": "type",
    "text": "Enumeration pd defining whether p or d are known for the boundary condition\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#Types-1",
    "page": "Types",
    "title": "Types",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nOrder   = [:type]"
},

{
    "location": "lib/Functions.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Functions.html#ModiaMedia.density-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.density",
    "category": "method",
    "text": "d = density(state) - return density from state in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_pTX",
    "category": "method",
    "text": "d = density_pTX(medium,p,T,X) - return density for medium from p, T, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_ph",
    "category": "method",
    "text": "d = density_ph(medium,p,h) - return density in [kg/m^3] for medium::PureSubstance from pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_phX",
    "category": "method",
    "text": "d = density_phX(medium,p,h,X) - return density for medium from p, h, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.dynamicViscosity-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.dynamicViscosity",
    "category": "method",
    "text": "eta = dynamicViscosity(state) - return dynamic viscosity state in [Pa*s]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.getMedium-Tuple{AbstractString}",
    "page": "Functions",
    "title": "ModiaMedia.getMedium",
    "category": "method",
    "text": "Medium = getMedium(name::AbstractString)\n\nReturn Medium object from medium name.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.isenthalpicState!-Tuple{ThermodynamicState,ThermodynamicState,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.isenthalpicState!",
    "category": "method",
    "text": "isenthalpicState!(state_b,state_a,dp)\n\nUpdate stateb by an isenthalpic transformation of statea with pressure drop dp.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.isenthalpicState-Tuple{AbstractMedium,ThermodynamicState,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.isenthalpicState",
    "category": "method",
    "text": "state_b = isenthalpicState(state_a,dp)\n\nReturn stateb by an isenthalpic transformation of statea with pressure drop dp.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.listMedia-Tuple{}",
    "page": "Functions",
    "title": "ModiaMedia.listMedia",
    "category": "method",
    "text": "listMedia()\n\nList available media of ModiaMedia\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.pressure",
    "category": "method",
    "text": "p = pressure(state) - return pressure from state in [Pa]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.pressure_dT",
    "category": "method",
    "text": "p = pressure_dT(medium,d,T) - return pressure in [Pa] for medium::PureSubstance from density c in [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dT!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dT!",
    "category": "method",
    "text": "setState_dT!(state, d,T)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  density d [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dT",
    "category": "method",
    "text": "state = setState_dT(medium, d,T)\n\nGenerate a state object for medium medium::PureSubstance for density d [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dTX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dTX!",
    "category": "method",
    "text": "setState_dTX!(state, d,T,X)\n\nUpdate the state::ThermodynamicState object with  density d [kg/m^3], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dTX",
    "category": "method",
    "text": "state = setState_dTX(medium, d,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for density d [kg/m^3], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pT!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pT!",
    "category": "method",
    "text": "setState_pT!(state, p,T)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pT-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pT",
    "category": "method",
    "text": "state = setState_pT(medium, p,T)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pTX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pTX!",
    "category": "method",
    "text": "setState_pTX!(state, p,T,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pTX",
    "category": "method",
    "text": "state = setState_pTX(medium, p,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ph!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ph!",
    "category": "method",
    "text": "setState_ph!(state, p,h)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and specific enthalpy h [J/kg]].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ph",
    "category": "method",
    "text": "state = setState_ph(medium, p,h)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_phX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_phX!",
    "category": "method",
    "text": "setState_phX!(state, p,h,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], specific enthalpy h [J/kg]] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_phX",
    "category": "method",
    "text": "state = setState_phX(medium, p,h,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific enthalpy h [J/kg]] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ps!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ps!",
    "category": "method",
    "text": "setState_ps!(state, p,s)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and specific entropy s [J/(kg*K)].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ps-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ps",
    "category": "method",
    "text": "state = setState_ps(medium, p,s)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific entropy s [J/(kg*K)].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_psX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_psX!",
    "category": "method",
    "text": "setState_psX!(state, p,s,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], specific entropy s [J/(kg*K)] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_psX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_psX",
    "category": "method",
    "text": "state = setState_psX(medium, p,s,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific entropy s [J/(kg*K)] and mass fractions vector X  or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy",
    "category": "method",
    "text": "h = specificEnthalpy(state) - return specific enthalpy from state in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy_dT",
    "category": "method",
    "text": "h = specificEnthalpy_dT(medium,d,T) - return specific enthalpy in [J/kg] for medium from density d [kg/m^3] and and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy_pTX",
    "category": "method",
    "text": "h = specificEnthalpy_pTX(medium,p,T,X) - return specific enthalpy for medium from p, T, and X or Xi in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificHeatCapacityCp-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificHeatCapacityCp",
    "category": "method",
    "text": "cp = specificHeatCapacityCp(state) - return specific heat capacity at constant pressure from state in [J/(kg*K)]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificInternalEnergy-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificInternalEnergy",
    "category": "method",
    "text": "u = specificInternalEnergy(state) - return specific internal energy at state in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.temperature",
    "category": "method",
    "text": "T = temperature(state) - return temperature from state in [K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.temperature_ph",
    "category": "method",
    "text": "T = temperature_ph(medium,p,h) - return temperature in [K] for medium::PureSubstance from pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.temperature_phX",
    "category": "method",
    "text": "T = temperature_phX(medium,p,h,X) - return temperature for medium from p, h, and X or Xi in [K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.cp_T-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.cp_T",
    "category": "method",
    "text": "cp = cp_T(data::SingleGasNasaData, T) - Compute specific heat capacity at constant pressure from temperature and gas data\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.dynamicViscosityLowPressure-NTuple{6,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.dynamicViscosityLowPressure",
    "category": "method",
    "text": "Dynamic viscosity of low pressure gases\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.enthalpyOfWater-Tuple{Float64}",
    "page": "Functions",
    "title": "ModiaMedia.enthalpyOfWater",
    "category": "method",
    "text": "Return specific enthalpy of water (solid/liquid) near atmospheric pressure from temperature T\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.h_T-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.h_T",
    "category": "method",
    "text": "h = h_T(data, T, exclEnthForm=true, refChoice=ReferenceEnthalpy_ZeroAt0K, h_off=0.0)\n\nReturn specific enthalpy from temperature and gas data.\n\nArguments\n\ndata::SingleGasNasaData: Data of the SingleGasNasa medium.\nT::Float64: Temperature in [K].\nexclEnthForm::Bool: If true, enthalpy of formation Hf is not included in specific enthalpy h.\nrefChoice::ReferenceEnthalpy: Choice of reference enthalpy.\nh_off::Float64: User defined offset for reference enthalpy, if SingleGasNasareferenceEnthalpy = ReferenceEnthalpyUserDefined\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.h_Tlow-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.h_Tlow",
    "category": "method",
    "text": "Compute specific enthalpy, low T region; reference is decided by the refChoice input, or by the referenceChoice package constant by default\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.s0_T-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.s0_T",
    "category": "method",
    "text": "Compute specific entropy from temperature and gas data\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.saturationPressure-Tuple{Float64}",
    "page": "Functions",
    "title": "ModiaMedia.saturationPressure",
    "category": "method",
    "text": "Return saturation pressure of water as a function of temperature T between 190 and 647.096 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.saturationPressureLiquid-Tuple{Float64}",
    "page": "Functions",
    "title": "ModiaMedia.saturationPressureLiquid",
    "category": "method",
    "text": "Return saturation pressure of water as a function of temperature T in the range of 273.16 to 647.096 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.spliceFunction-NTuple{4,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.spliceFunction",
    "category": "method",
    "text": "Spline interpolation of two functions\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardCharacteristics-Tuple{AbstractMedium}",
    "page": "Functions",
    "title": "ModiaMedia.standardCharacteristics",
    "category": "method",
    "text": "dict = standardCharacteristics(medium::AbstractMedium)\n\nReturn a dict::dict{AbstractString,Any} dictionary with the most  important characteristics of the medium as vectors.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardPlot-Tuple{AbstractMedium}",
    "page": "Functions",
    "title": "ModiaMedia.standardPlot",
    "category": "method",
    "text": "standardPlot(medium::AbstractMedium; figure=1)\n\nPlot the standardCharacteristics(medium) of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.sublimationPressureIce-Tuple{Float64}",
    "page": "Functions",
    "title": "ModiaMedia.sublimationPressureIce",
    "category": "method",
    "text": "Return sublimation pressure of water as a function of temperature T between 190 and 273.16 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nOrder   = [:function]"
},

]}
