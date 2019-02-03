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
    "text": "  using ModiaMedia\r\n\r\n  # Define medium to be used\r\n  Medium = getMedium(\"N2\");\r\n\r\n  # Define the operating point where the medium shall be evaluated.\r\n  p = 1e5    # in [Pa]\r\n  T = 300.0  # in [K]\r\n\r\n  # Set the medium-specific thermodynamic state from p and T\r\n  # (could be also set from p,h, or p,s, or d,T, or\r\n  # p,T,X, or p,h,X, or p,s,X, or d,T,X)\r\n  state = setState_pT(Medium, p, T)\r\n\r\n  # Update a state object with new values\r\n  setState_pT!(state, p, T)\r\n\r\n  # Call media functions (here to compute density and specific enthalpy)\r\n  d = density(state)\r\n  h = specificEnthalpy(state)\r\n\r\n  # Print computed values\r\n  println(\"data for p=$p, T=$T:\")\r\n  println(\"density          = \", d)\r\n  println(\"specificEnthalpy = \", h)\r\n\r\n  # List the available media\r\n  listMedia()\r\n\r\n  # Plot the most important characteristics of the medium\r\n  standardPlot(Medium)This example generates the following plot:(Image: standardPlot)"
},

{
    "location": "index.html#Available-media-1",
    "page": "Home",
    "title": "Available media",
    "category": "section",
    "text": "The available media can be listed with listMedia() resulting in:Row name type\n1 Ar SingleGasNasa\n2 C2H2_vinylidene SingleGasNasa\n3 C2H4 SingleGasNasa\n4 C2H5OH SingleGasNasa\n5 C2H6 SingleGasNasa\n6 C3H6_propylene SingleGasNasa\n7 C3H8 SingleGasNasa\n8 C4H10_n_butane SingleGasNasa\n9 C4H8_1_butene SingleGasNasa\n10 C5H10_1_pentene SingleGasNasa\n11 C5H12_n_pentane SingleGasNasa\n12 C6H12_1_hexene SingleGasNasa\n13 C6H14_n_hexane SingleGasNasa\n14 C6H6 SingleGasNasa\n15 C7H14_1_heptene SingleGasNasa\n16 C7H16_n_heptane SingleGasNasa\n17 C8H10_ethylbenz SingleGasNasa\n18 C8H18_n_octane SingleGasNasa\n19 CH3OH SingleGasNasa\n20 CH4 SingleGasNasa\n21 CL2 SingleGasNasa\n22 CO SingleGasNasa\n23 CO2 SingleGasNasa\n24 ConstantPropertyLiquidWater SimpleMedium\n25 F2 SingleGasNasa\n26 H2 SingleGasNasa\n27 H2O SingleGasNasa\n28 He SingleGasNasa\n29 MoistAir MoistAir\n30 N2 SingleGasNasa\n31 N2O SingleGasNasa\n32 NH3 SingleGasNasa\n33 NO SingleGasNasa\n34 NO2 SingleGasNasa\n35 Ne SingleGasNasa\n36 O2 SingleGasNasa\n37 SO2 SingleGasNasa\n38 SO3 SingleGasNasa\n39 SimpleAir SimpleIdealGasMedium"
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
    "text": "Martin Otter (DLR - Institute of System Dynamics and Control)\n\nHilding Elmqvist (Mogram),\n\nChris Laughman (MERL).\nAll the content of ModiaMedia is based on Modelica.Media which was and is developed from many people.License: MIT (expat)"
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
    "location": "lib/MediaTypes.html#",
    "page": "Media Types",
    "title": "Media Types",
    "category": "page",
    "text": ""
},

{
    "location": "lib/MediaTypes.html#Media-Types-1",
    "page": "Media Types",
    "title": "Media Types",
    "category": "section",
    "text": "In this section the currently available media types are summarized. The Available media list the actually available media. A concrete medium object is inquired with getMedium(mediumName)."
},

{
    "location": "lib/MediaTypes.html#SimpleMedium-1",
    "page": "Media Types",
    "title": "SimpleMedium",
    "category": "section",
    "text": "PureSubstance Medium model with linear dependency of specific internal energy and specific enthalpy from temperature. All other quantities, especially density, are constant. Example:# Definition of medium in ModiaMedia/dict/SimpleMedium.jl\r\ndict[\"ConstantPropertyLiquidWater\"] = ModiaMedia.SimpleMedium(\r\n                        mediumName = \"ConstantPropertyLiquidWater\",\r\n\r\n                        fluidConstants = ModiaMedia.BasicFluidConstants(\r\n                            chemicalFormula=\"H2O\", \r\n                            structureFormula=\"H2O\", \r\n                            casRegistryNumber=\"7732-18-5\", \r\n                            iupacName=\"oxidane\", \r\n                            molarMass=0.018015268),\r\n\r\n                        data = ModiaMedia.SimpleMediumData(\r\n                            cp_const=4184, \r\n                            cv_const=4184, \r\n                            d_const=995.586, \r\n                            eta_const=1.e-3, \r\n                            lambda_const=0.598, \r\n                            a_const=1484, \r\n                            T_min=ModiaMedia.from_degC(-1), \r\n                            T_max=ModiaMedia.from_degC(130), \r\n                            T0=273.15, \r\n                            MM_const=0.018015268)\r\n                    )\r\n\r\n\r\n# Using this medium\r\nusing ModiaMedium\r\nstandardPlot( getMedium(\"ConstantPropertyLiquidWater\") )(Image: standardPlot)"
},

{
    "location": "lib/MediaTypes.html#SimpleIdealGasMedium-1",
    "page": "Media Types",
    "title": "SimpleIdealGasMedium",
    "category": "section",
    "text": "PureSubstance Medium model of ideal gas with linear dependency of specific  internal energy and specific enthalpy from temperature and constant transport  properties. Density is a function of temperature and pressure. Example:# Definition of medium in ModiaMedia/dict/SimpleIdealGasMedium.jl\r\ndict[\"SimpleAir\"] = ModiaMedia.SimpleIdealGasMedium(\r\n                        mediumName = \"SimpleAir\",\r\n\r\n                        fluidConstants = ModiaMedia.BasicFluidConstants(\r\n                            iupacName=\"simple air\",\r\n                            casRegistryNumber=\"not a real substance\",\r\n                            chemicalFormula=\"N2, O2\",\r\n                            structureFormula=\"N2, O2\",\r\n                            molarMass=0.0289651159),\r\n\r\n                        data = ModiaMedia.SimpleIdealGasMediumData(\r\n                            cp_const     = 1005.45, \r\n                            MM_const     = 0.0289651159, \r\n                            R_gas        = 8.3144598/0.0289651159, \r\n                            eta_const    = 1.82e-5, \r\n                            lambda_const = 0.026, \r\n                            T_min        = ModiaMedia.from_degC(0.0), \r\n                            T_max        = ModiaMedia.from_degC(100.0),\r\n                            T0           = 298.15)\r\n                    )\r\n\r\n\r\n# Using this medium\r\nusing ModiaMedium\r\nstandardPlot( getMedium(\"SimpleAir\") )(Image: standardPlot)"
},

{
    "location": "lib/MediaTypes.html#SingleGasNasa-1",
    "page": "Media Types",
    "title": "SingleGasNasa",
    "category": "section",
    "text": "PureSubstance medium model of ideal gas with nonlinear functions for specific enthalpy and specific internal energy.  Density is a function of temperature and pressure, whereas all other functions are a function of temperature only. The nonlinear function for specific enthalpy and its data are fromMcBride B.J., Zehe M.J., and Gordon S. (2002):  NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species.  NASA report TP-2002-211556# Using this medium\r\nusing ModiaMedium\r\nstandardPlot( getMedium(\"N2\") )(Image: standardPlot)"
},

{
    "location": "lib/MediaTypes.html#MoistAir-1",
    "page": "Media Types",
    "title": "MoistAir",
    "category": "section",
    "text": ""
},

{
    "location": "lib/MediaTypes.html#Thermodynamic-Model-1",
    "page": "Media Types",
    "title": "Thermodynamic Model",
    "category": "section",
    "text": "CondensingGases thermodynamic model of moist air including the  fog region and temperatures below zero degC. The governing assumptions in this model are:the perfect gas law applies\nwater volume other than that of steam is neglectedAll extensive properties are expressed in terms of the total mass in order to comply with other media in this library. However, for moist air it is rather common to express the absolute humidity in terms of mass of dry air only, which has advantages when working with charts. In addition, care must be taken, when working with mass fractions with respect to total mass, that all properties refer to the same water content when being used in mathematical operations (which is always the case if based on dry air only). Therefore two absolute humidities are computed in the BaseProperties model: X denotes the absolute humidity in terms of the total mass while x denotes the absolute humidity per unit mass of dry air. In addition, the relative humidity phi is also computed.At the triple point temperature of water of 0.01 °C or 273.16 K and a relative humidity greater than 1 fog may be present as liquid and as ice resulting in a specific enthalpy somewhere between those of the two isotherms for solid and liquid fog, respectively. For numerical reasons a coexisting mixture of 50% solid and 50% liquid fog is assumed in the fog region at the triple point in this model. "
},

{
    "location": "lib/MediaTypes.html#Range-of-validity-1",
    "page": "Media Types",
    "title": "Range of validity",
    "category": "section",
    "text": "From the assumptions mentioned above it follows that the pressure should be in the region around atmospheric conditions or below (a few bars may still be fine though). Additionally a very high water content at low temperatures would yield incorrect densities, because the volume of the liquid or solid phase would not be negligible anymore. The model does not provide information on limits for water drop size in the fog region or transport information for the actual condensation or evaporation process in combination with surfaces. All excess water which is not in its vapour state is assumed to be still present in the air regarding its energy but not in terms of its spatial extent.The thermodynamic model may be used for temperatures ranging from 190 ... 647 K. This holds for all functions unless otherwise stated in their description. However, although the model works at temperatures above the saturation temperature it is questionable to use the term \"relative humidity\" in this region. Please note, that although several functions compute pure water properties, they are designed to be used within the moist air medium model where properties are dominated by air and steam in their vapor states, and not for pure liquid water applications. "
},

{
    "location": "lib/MediaTypes.html#Transport-Properties-1",
    "page": "Media Types",
    "title": "Transport Properties",
    "category": "section",
    "text": "Several additional functions that are not needed to describe the thermodynamic system, but are required to model transport processes, like heat and mass transfer, may be called. They usually neglect the moisture influence unless otherwise stated. "
},

{
    "location": "lib/MediaTypes.html#Application-1",
    "page": "Media Types",
    "title": "Application",
    "category": "section",
    "text": "The model\'s main area of application is all processes that involve moist air cooling under near atmospheric pressure with possible moisture condensation. This is the case in all domestic and industrial air conditioning applications. Another large domain of moist air applications covers all processes that deal with dehydration of bulk material using air as a transport medium. "
},

{
    "location": "lib/MediaTypes.html#Usage-1",
    "page": "Media Types",
    "title": "Usage",
    "category": "section",
    "text": "# Using this medium\r\nusing ModiaMedium\r\nstandardPlot( getMedium(\"MoistAir\") )(Image: standardPlot)"
},

{
    "location": "lib/Types.html#",
    "page": "Exported Types",
    "title": "Exported Types",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Types.html#Exported-Types-1",
    "page": "Exported Types",
    "title": "Exported Types",
    "category": "section",
    "text": ""
},

{
    "location": "lib/Types.html#Index-1",
    "page": "Exported Types",
    "title": "Index",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nPrivate = false\r\nOrder   = [:type]"
},

{
    "location": "lib/Types.html#ModiaMedia.AbstractFluidConstants",
    "page": "Exported Types",
    "title": "ModiaMedia.AbstractFluidConstants",
    "category": "type",
    "text": "abstract type AbstractFluidConstants\n\nAbstract type of all FluidConstants structures.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.AbstractMedium",
    "page": "Exported Types",
    "title": "ModiaMedia.AbstractMedium",
    "category": "type",
    "text": "abstract type AbstractMedium\n\nAbstract type of all media.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.BasicFluidConstants",
    "page": "Exported Types",
    "title": "ModiaMedia.BasicFluidConstants",
    "category": "type",
    "text": "fluidConstants = BasicFluidConstants(;iupacName=\"\", casRegistryNumber=\"\", \n                     chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN)\n\nGenerate a BasicFluidConstants <: AbstractFluidConstants object  containing the minimal information about the standard data of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.CondensingGases",
    "page": "Exported Types",
    "title": "ModiaMedia.CondensingGases",
    "category": "type",
    "text": "abstract type CondensingGases <: AbstractMedium\n\nAbstract type of all media consisting of condensing media.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidInfos",
    "page": "Exported Types",
    "title": "ModiaMedia.FluidInfos",
    "category": "type",
    "text": "infos = FluidInfos(; <keyword arguments, see below>)\n\nGenerate a new FluidInfos object, containing generic properties of a medium.\n\nKeyword arguments\n\nName and type Description\nmediumName::AbstractString Name of the medium\nsubstanceNames::Vector{AbstractString} Names of the mixture substances. Set substanceNames=[mediumName] if only one substance.\nextraPropertiesNames::Vector{AbstractString} Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused\nThermoStates::IndependentVariables Enumeration type for independent variables\nsingleState::Bool = true, if u and d are not a function of pressure\nreducedX::Bool = true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance\nfixedX::Bool = true if medium contains the equation X = reference_X\nreference_p Reference pressure of Medium: default 1 atmosphere\nreference_T Reference temperature of Medium: default 25 deg Celsius\nreference_X Reference mass fractions of medium\np_default Default value for pressure of medium (for initialization)\nT_default Default value for temperature of medium (for initialization)\nh_default Default value for specific enthalpy of medium (for initialization)\nX_default Default value for specific enthalpy of medium (for initialization)\nnS Number of substances\nnX Number of mass fractions\nnXi Number of structurally independent mass fractions\nnC Number of extra (outside of standard mass-balance) transported properties\nC_nominal Default for the nominal values for the extra properties\n\nExample\n\nimport ModiaMedia\n\ninfos = ModiaMedia.FluidInfos(mediumName           = \"simpleMedium\",\n                              substanceNames       = [mediumName],\n                              extraPropertiesNames = fill(\"\",0),\n                              ThermoStates         = IndependentVariables_T)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidLimits",
    "page": "Exported Types",
    "title": "ModiaMedia.FluidLimits",
    "category": "type",
    "text": "fluidLimits = FluidLimits(; TMIN=NaN, TMAX=NaN, DMIN=NaN, DMAX=NaN, PMIN=NaN, PMAX=NaN, \n                            HMIN=NaN, HMAX=NaN, SMIN=NaN, SMAX=NaN)\n\nGenerate a new FluidLimits object, containing the validity limits of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IdealGasFluidConstants",
    "page": "Exported Types",
    "title": "ModiaMedia.IdealGasFluidConstants",
    "category": "type",
    "text": "fluidConstants = IdealGasFluidConstants(; iupacName=\"\", casRegistryNumber=\"\", \n      chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN,\n      criticalTemperature=NaN, criticalPressure=NaN, criticalMolarVolume=NaN, \n      acentricFactor=NaN, meltingPoint=NaN, normalBoilingPoint=NaN, dipoleMoment=NaN,\n      hasIdealGasHeatCapacity=false, hasCriticalData=false, hasDipoleMoment=false,\n      hasFundamentalEquation=false, hasLiquidHeatCapacity=false, hasSolidHeatCapacity=false,\n      hasAccurateViscosityData=false, hasAccurateConductivityData=false,\n      hasVapourPressureCurve=false, hasAcentricFactor=false,\n      HCRIT0=0.0, SCRIT0=0.0, deltah=0.0, deltas=0.0)\n\nGenerate a IdealGasFluidConstants <: AbstractFluidConstants object containing the minimal information about the standard data of ideal gas media (critical, triple, molecular and other standard data).\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IndependentVariables",
    "page": "Exported Types",
    "title": "ModiaMedia.IndependentVariables",
    "category": "type",
    "text": "@enum IndependentVariables\n\nEnumeration defining the independent variables of a medium. Possible values:\n\nValue Independent variables\nIndependentVariables_T Temperature\nIndependentVariables_pT Pressure, temperature\nIndependentVariables_ph Pressure, specific enthalpy\nIndependentVariables_phX Pressure, specific enthalpy, mass fractions\nIndependentVariables_pTX Pressure, temperature, mass fractions\nIndependentVariables_dTX Density, temperature, mass fractions\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MixtureMedium",
    "page": "Exported Types",
    "title": "ModiaMedia.MixtureMedium",
    "category": "type",
    "text": "abstract type MixtureMedium <: AbstractMedium\n\nAbstract type of all media consisting of a mixture of media.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MixtureThermodynamicState",
    "page": "Exported Types",
    "title": "ModiaMedia.MixtureThermodynamicState",
    "category": "type",
    "text": "abstract type MixtureThermodynamicState <: ThermodynamicState\n\nAbstract type of the states of all media consisting of a mixture of media.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.MoistAirState",
    "page": "Exported Types",
    "title": "ModiaMedia.MoistAirState",
    "category": "type",
    "text": "state = MoistAirState(Medium, p, T, X)\n\nGenerate a MoistAirState <: MixtureThermodynamicState object containing pressure p [Pa], temperature T [K], and a vector of mass fractions as states, where X[1] is the mass fraction of Steam and X[2] is the mass fraction of dry air. If argument X has only one element, X[2] is computed from X[1] and stored in the state.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.PureSubstance",
    "page": "Exported Types",
    "title": "ModiaMedia.PureSubstance",
    "category": "type",
    "text": "abstract type PureSubstance <: AbstractMedium\n\nAbstract type of all media consisting of a pure substance.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ReferenceEnthalpy",
    "page": "Exported Types",
    "title": "ModiaMedia.ReferenceEnthalpy",
    "category": "type",
    "text": "@enum ReferenceEnthalpy\n\nEnumeration defining the reference enthalpy of a medium. Possible values:\n\nReferenceEnthalpy_ZeroAt0K: The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded\nReferenceEnthalpy_ZeroAt25C: The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded\nReferenceEnthalpy_UserDefined: The user-defined reference enthalpy is used at 293.15 K (25 degC)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ReferenceEntropy",
    "page": "Exported Types",
    "title": "ModiaMedia.ReferenceEntropy",
    "category": "type",
    "text": "@enum ReferenceEntropy\n\nEnumeration defining the reference entropy of a medium. Possible values:\n\nReferenceEntropy_ZeroAt0K: The entropy is 0 at 0 K (default)\nReferenceEntropy_ZeroAt0C: The entropy is 0 at 0 degC\nReferenceEntropy_UserDefined: The user-defined reference entropy is used at 293.15 K (25 degC)\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleIdealGasMediumState",
    "page": "Exported Types",
    "title": "ModiaMedia.SimpleIdealGasMediumState",
    "category": "type",
    "text": "state = SimpleIdealGasMediumState(Medium, p, T)\n\nGenerate a SimpleIdealGasMediumState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleMediumState",
    "page": "Exported Types",
    "title": "ModiaMedia.SimpleMediumState",
    "category": "type",
    "text": "state = SimpleMediumState(Medium, p, T)\n\nGenerate a SimpleMediumState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasaState",
    "page": "Exported Types",
    "title": "ModiaMedia.SingleGasNasaState",
    "category": "type",
    "text": "state = SingleGasNasaState(Medium, p, T)\n\nGenerate a SingleGasNasaState <: ThermodynamicState object containing pressure p [Pa] and temperature T [K] as thermodynamic states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.ThermodynamicState",
    "page": "Exported Types",
    "title": "ModiaMedia.ThermodynamicState",
    "category": "type",
    "text": "abstract type ThermodynamicState\n\nAbstract type of all media states.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#Documentation-1",
    "page": "Exported Types",
    "title": "Documentation",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nPrivate = false\r\nOrder   = [:type]"
},

{
    "location": "lib/Functions.html#",
    "page": "Exported Functions",
    "title": "Exported Functions",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Functions.html#Exported-Functions-1",
    "page": "Exported Functions",
    "title": "Exported Functions",
    "category": "section",
    "text": ""
},

{
    "location": "lib/Functions.html#Index-1",
    "page": "Exported Functions",
    "title": "Index",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nPrivate = false\r\nOrder   = [:function]"
},

{
    "location": "lib/Functions.html#ModiaMedia.density-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.density",
    "category": "method",
    "text": "density(state)\n\nReturn density from state::ThermodynamicState in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.density_pTX",
    "category": "method",
    "text": "density_pTX(medium,p,T,X)\n\nReturn density for medium::AbstractMedium from p, T, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.density_ph",
    "category": "method",
    "text": "density_ph(medium,p,h)\n\nReturn density in [kg/m^3] for medium::PureSubstance from  pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.density_phX",
    "category": "method",
    "text": "density_phX(medium,p,h,X)\n\nReturn density for medium::AbstractMedium from p, h, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.dynamicViscosity-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.dynamicViscosity",
    "category": "method",
    "text": "dynamicViscosity(state)\n\nReturn dynamic viscosity from state::ThermodynamicState in [Pa*s]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.gasConstant-Tuple{MixtureMedium,MixtureThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.gasConstant",
    "category": "method",
    "text": "R = gasConstant(state)\n\nReturn gas constant of MixtureMedium from state::MixtureThermodynamicState in [J/mol.K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.getMedium-Tuple{AbstractString}",
    "page": "Exported Functions",
    "title": "ModiaMedia.getMedium",
    "category": "method",
    "text": "Medium = getMedium(name::AbstractString)\n\nReturn Medium object from medium name.  Possible values of argument name can be inquired via listMedia() (Available media).\n\nExamples\n\nmedium = getMedium(\"SimpleAir\")\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.isenthalpicState!-Tuple{ThermodynamicState,ThermodynamicState,Float64}",
    "page": "Exported Functions",
    "title": "ModiaMedia.isenthalpicState!",
    "category": "method",
    "text": "isenthalpicState!(state_b,state_a,dp)\n\nUpdate state_b by an isenthalpic transformation of state_a with pressure drop dp:\n\npressure(state_b)         = pressure(state_a) + dp\nspecificEnthalpy(state_b) = specificEnthalpy(state_a)\nstate_b.X                 = state_a.X\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.isenthalpicState-Tuple{AbstractMedium,ThermodynamicState,Float64}",
    "page": "Exported Functions",
    "title": "ModiaMedia.isenthalpicState",
    "category": "method",
    "text": "state_b = isenthalpicState(state_a,dp)\n\nReturn state_b by an isenthalpic transformation of state_a with pressure drop dp:\n\npressure(state_b)         = pressure(state_a) + dp\nspecificEnthalpy(state_b) = specificEnthalpy(state_a)\nstate_b.X                 = state_a.X\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.listMedia-Tuple{}",
    "page": "Exported Functions",
    "title": "ModiaMedia.listMedia",
    "category": "method",
    "text": "listMedia()\n\nList available media of ModiaMedia.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.pressure",
    "category": "method",
    "text": "pressure(state)\n\nReturn pressure from state::ThermodynamicState in [Pa]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.pressure_dT",
    "category": "method",
    "text": "pressure_dT(medium,d,T)\n\nReturn pressure in [Pa] for medium::PureSubstance from density d in [kg/m^3]  and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.saturationPressureOfLiquidWater-Tuple{Float64}",
    "page": "Exported Functions",
    "title": "ModiaMedia.saturationPressureOfLiquidWater",
    "category": "method",
    "text": "saturationPressureOfLiquidWater(Tsat)\n\nReturn saturation pressure of liquid water as a function of saturation temperature Tsat in the range of 273.16 to 647.096 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.saturationPressureOfWater-Tuple{Float64}",
    "page": "Exported Functions",
    "title": "ModiaMedia.saturationPressureOfWater",
    "category": "method",
    "text": "saturationPressureOfWater(Tsat)\n\nReturn saturation pressure of water as a function of  saturation temperature Tsat between 190 K and 647.096 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dT!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_dT!",
    "category": "method",
    "text": "setState_dT!(state, d,T)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  density d [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_dT",
    "category": "method",
    "text": "state = setState_dT(medium, d,T)\n\nGenerate a state object for medium medium::PureSubstance for density d [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dTX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_dTX!",
    "category": "method",
    "text": "setState_dTX!(state, d,T,X)\n\nUpdate the state::ThermodynamicState object with  density d [kg/m^3], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_dTX",
    "category": "method",
    "text": "state = setState_dTX(medium, d,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for density d [kg/m^3], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pT!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_pT!",
    "category": "method",
    "text": "setState_pT!(state, p,T)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pT-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_pT",
    "category": "method",
    "text": "state = setState_pT(medium, p,T)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pTX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_pTX!",
    "category": "method",
    "text": "setState_pTX!(state, p,T,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_pTX",
    "category": "method",
    "text": "state = setState_pTX(medium, p,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ph!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_ph!",
    "category": "method",
    "text": "setState_ph!(state, p,h)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and specific enthalpy h [J/kg]].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_ph",
    "category": "method",
    "text": "state = setState_ph(medium, p,h)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_phX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_phX!",
    "category": "method",
    "text": "setState_phX!(state, p,h,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], specific enthalpy h [J/kg]] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_phX",
    "category": "method",
    "text": "state = setState_phX(medium, p,h,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific enthalpy h [J/kg]] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ps!-Tuple{ThermodynamicState,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_ps!",
    "category": "method",
    "text": "setState_ps!(state, p,s)\n\nUpdate the state::ThermodynamicState of a PureSubstance medium with  pressure p [Pa] and specific entropy s [J/(kg*K)].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ps-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_ps",
    "category": "method",
    "text": "state = setState_ps(medium, p,s)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific entropy s [J/(kg*K)].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_psX!-Tuple{ThermodynamicState,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_psX!",
    "category": "method",
    "text": "setState_psX!(state, p,s,X)\n\nUpdate the state::ThermodynamicState object with  pressure p [Pa], specific entropy s [J/(kg*K)] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_psX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.setState_psX",
    "category": "method",
    "text": "state = setState_psX(medium, p,s,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific entropy s [J/(kg*K)] and mass fractions vector X  or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.specificEnthalpy",
    "category": "method",
    "text": "specificEnthalpy(state)\n\nReturn specific enthalpy from state::ThermodynamicState in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_dT-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.specificEnthalpy_dT",
    "category": "method",
    "text": "specificEnthalpy_dT(medium,d,T)\n\nReturn specific enthalpy in [J/kg] for medium::PureSubstance from density d [kg/m^3] and and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_pTX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.specificEnthalpy_pTX",
    "category": "method",
    "text": "specificEnthalpy_pTX(medium,p,T,X)\n\nReturn specific enthalpy for medium::AbstractMedium from p, T, and X or Xi in [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificHeatCapacityCp-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.specificHeatCapacityCp",
    "category": "method",
    "text": "specificHeatCapacityCp(state)\n\nReturn specific heat capacity at constant pressure from state::ThermodynamicState in [J/(kg*K)]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificInternalEnergy-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.specificInternalEnergy",
    "category": "method",
    "text": "specificInternalEnergy(state)\n\nReturn specific internal energy from state::ThermodynamicState in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardCharacteristics-Tuple{AbstractMedium}",
    "page": "Exported Functions",
    "title": "ModiaMedia.standardCharacteristics",
    "category": "method",
    "text": "dict = standardCharacteristics(medium::AbstractMedium)\n\nReturn a dict::dict{AbstractString,Any} dictionary with the most  important characteristics of the medium as vectors.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardPlot-Tuple{AbstractMedium}",
    "page": "Exported Functions",
    "title": "ModiaMedia.standardPlot",
    "category": "method",
    "text": "standardPlot(medium::AbstractMedium; figure=1)\n\nPlot the standardCharacteristics(medium) of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.sublimationPressureIce-Tuple{Float64}",
    "page": "Exported Functions",
    "title": "ModiaMedia.sublimationPressureIce",
    "category": "method",
    "text": "sublimationPressureIce(Tsat)\n\nReturn sublimation pressure of water as a function of saturation temperature Tsat  between 190 and 273.16 K\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature-Tuple{AbstractMedium,ThermodynamicState}",
    "page": "Exported Functions",
    "title": "ModiaMedia.temperature",
    "category": "method",
    "text": "temperature(state)\n\nReturn temperature from state::ThermodynamicState in [K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_ph-Tuple{PureSubstance,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.temperature_ph",
    "category": "method",
    "text": "temperature_ph(medium,p,h)\n\nReturn temperature in [K] for medium::PureSubstance from  pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_phX-Tuple{AbstractMedium,Any,Any,Any}",
    "page": "Exported Functions",
    "title": "ModiaMedia.temperature_phX",
    "category": "method",
    "text": "temperature_phX(medium,p,h,X)\n\nReturn temperature for medium::AbstractMedium from p, h, and X or Xi in [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#Documentation-1",
    "page": "Exported Functions",
    "title": "Documentation",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nPrivate = false\r\nOrder   = [:function]"
},

]}
