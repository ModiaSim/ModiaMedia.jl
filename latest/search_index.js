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
    "text": "ModiaMedia shall provide Media models  for use with Modia and other Julia packages. The initial goal is to achieve a similar functionality as Modelica.Media, the standard media library for Modelica models, but with improvements based on Julia features such as multiple dispatch.This package is under development and it is planned to provide all media from Modelica.Media in this package."
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
    "text": "  import ModiaMedia\r\n\r\n  # Define medium to be used\r\n  medium = ModiaMedia.Medium(\"N2\");\r\n\r\n  # Define the operating point where the medium shall be evaluated. \r\n  p = 1e5    # in [Pa]\r\n  T = 300.0  # in [K]\r\n\r\n  # Determine the medium-specific thermodynamic state\r\n  # (setState_ph, setState_ps, setState_dT could also be used)\r\n  state = ModiaMedia.setState_pT(medium, p, T)\r\n\r\n  # Call media functions (here to compute density and specific enthalpy)\r\n  d = ModiaMedia.density(medium,state)\r\n  h = ModiaMedia.specificEnthalpy(medium,state)\r\n\r\n  # Print computed values\r\n  println(\"data for p=$p, T=$T:\")\r\n  println(\"density          = \", d)\r\n  println(\"specificEnthalpy = \", h)\r\n\r\n  # Plot the most important characteristics of the medium\r\n  ModiaMedia.standardPlot(medium)This example generates the following plot:(Image: standardPlot)"
},

{
    "location": "index.html#Currently-available-media-1",
    "page": "Home",
    "title": "Currently available media",
    "category": "section",
    "text": "SimpleLiquidWaterThe following 37 ideal gases (from NASA Glenn coefficients):  Ar, CH4, CH3OH, CO, CO2, C2H2vinylidene, C2H4, C2H5OH, C2H6, C3H6propylene,  C3H8, C3H8O1propanol, C4H81butene, C4H10nbutane, C5H101pentene,  C5H12npentane, C6H6, C6H121hexene, C6H14nhexane, C7H141heptene,   C7H16nheptane, C8H10ethylbenz, C8H18noctane, CL2, F2, H2, H2O,     He, NH3, NO, NO2, N2, N2O, Ne, O2, SO2, SO3 "
},

{
    "location": "index.html#Structure-of-package-1",
    "page": "Home",
    "title": "Structure of package",
    "category": "section",
    "text": "A medium is a struct of the following type:struct MediumXXX <: AbstractMedium  # or of a subtype of AbstractMedium\r\n    infos::FluidInfos\r\n    constants::AbstractFluidConstants\r\n    limits::FluidLimits\r\n    data  # medium specific data\r\nendand all instances of this struct are stored in a dictionary. This dictionary is constructed in a preprocessing step by running \"ModiaMedia/dict/GenerateMediumDict.jl\". This module contains code that was mostly automatically converted from Modelica.Media to Julia. The resulting dictionary is serialized and stored in \"ModiaMedia/src/Media/media.julia_serializer\". When package ModiaMedia is compiled, this serialized dictionary is deserialized and included in the compiled package.Function ModiaMedia.Medium(name) returns the MediumXXX instance stored in the medium dictionary with key name."
},

{
    "location": "index.html#Status-1",
    "page": "Home",
    "title": "Status",
    "category": "section",
    "text": "The ModiaMedia package development has just started and a lot has to be improved."
},

{
    "location": "index.html#Release-Notes-1",
    "page": "Home",
    "title": "Release Notes",
    "category": "section",
    "text": ""
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
    "location": "lib/Types.html#ModiaMedia.AbstractFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.AbstractFluidConstants",
    "category": "type",
    "text": "abstract type AbstractFluidConstants - Abstract type of all FluidConstants structures\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.AbstractMedium",
    "page": "Types",
    "title": "ModiaMedia.AbstractMedium",
    "category": "type",
    "text": "abstract type AbstractMedium - Abstract type of all media\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.BasicFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.BasicFluidConstants",
    "category": "type",
    "text": "fluidConstants = BasicFluidConstants(;iupacName=\"\", casRegistryNumber=\"\", \n                     chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN)\n\nGenerate a BasicFluidConstants <: AbstractFluidConstants object  containing the minimal information about the standard data of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidInfos",
    "page": "Types",
    "title": "ModiaMedia.FluidInfos",
    "category": "type",
    "text": "infos = FluidInfos(;mediumName=Missing, substanceNames=[mediumName],\n                    extraPropertiesNames=fill(\"\",0), ThermoStates=Missing,\n                    singleState=Missing, reducedX=true, fixedX=false, \n                    reference_p=101325, reference_T=298.15,\n                    reference_X=fill(1/length(substanceNames),length(substanceNames)),\n                    p_default=101325, T_default=293.15, h_default=NaN,\n                    X_default=reference_X, C_nominal=1e-6*ones(length(extraPropertiesNames)))\n\nGenerate a new FluidInfos object defining the substance names, states, dimensions,  reference and default values for a medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.FluidLimits",
    "page": "Types",
    "title": "ModiaMedia.FluidLimits",
    "category": "type",
    "text": "limits = FluidLimits(; TMIN=NaN, TMAX=NaN, DMIN=NaN, DMAX=NaN, PMIN=NaN, PMAX=NaN, \n                       HMIN=NaN, HMAX=NaN, SMIN=NaN, SMAX=NaN)\n\nGenerate a new FluidLimits object, containing the validity limits of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IdealGasFluidConstants",
    "page": "Types",
    "title": "ModiaMedia.IdealGasFluidConstants",
    "category": "type",
    "text": "fluidConstants = IdealGasFluidConstants(; iupacName=\"\", casRegistryNumber=\"\", \n      chemicalFormula=\"\", structureFormula=\"\", molarMass=NaN,\n      criticalTemperature=NaN, criticalPressure=NaN, criticalMolarVolume=NaN, \n      acentricFactor=NaN, meltingPoint=NaN, normalBoilingPoint=NaN, dipoleMoment=NaN,\n      hasIdealGasHeatCapacity=false, hasCriticalData=false, hasDipoleMoment=false,\n      hasFundamentalEquation=false, hasLiquidHeatCapacity=false, hasSolidHeatCapacity=false,\n      hasAccurateViscosityData=false, hasAccurateConductivityData=false,\n      hasVapourPressureCurve=false, hasAcentricFactor=false,\n      HCRIT0=0.0, SCRIT0=0.0, deltah=0.0, deltas=0.0)\n\nGenerate a IdealGasFluidConstants <: AbstractFluidConstants object containing the minimal information about the standard data of ideal gas media (critical, triple, molecular and other standard data).\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.IndependentVariables",
    "page": "Types",
    "title": "ModiaMedia.IndependentVariables",
    "category": "type",
    "text": "@enum IndependentVariables\n\nEnumeration defining the independent variables of a medium. Possible values:\n\nValue Independent variables\nIndependentVariables_T Temperature\nIndependentVariables_pT Pressure, temperature\nIndependentVariables_ph Pressure, specific enthalpy\nIndependentVariables_phX Pressure, specific enthalpy, mass fractions\nIndependentVariables_pTX Pressure, temperature, mass fractions\nIndependentVariables_dTX Density, temperature, mass fractions\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.Init",
    "page": "Types",
    "title": "ModiaMedia.Init",
    "category": "type",
    "text": "Enumeration Init defining initialization for fluid flow\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.PureSubstance",
    "page": "Types",
    "title": "ModiaMedia.PureSubstance",
    "category": "type",
    "text": "abstract type PureSubstance <: AbstractMedium - Abstract type of all media consisting of a pure substance\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.PureSubstanceThermodynamicState",
    "page": "Types",
    "title": "ModiaMedia.PureSubstanceThermodynamicState",
    "category": "type",
    "text": "abstract type PureSubstanceThermodynamicState <: ThermodynamicState - Abstract type of the states of all media consisting of a pure substance\n\n\n\n\n\n"
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
    "location": "lib/Types.html#ModiaMedia.SimpleMedium",
    "page": "Types",
    "title": "ModiaMedia.SimpleMedium",
    "category": "type",
    "text": "medium = SimpleMedium(;infos=nothing, constants=nothing, limits=FluidLimits(), data=nothing)\n\nGenerate a SimpleMedium <: PureSubstance medium object.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SimpleMediumData",
    "page": "Types",
    "title": "ModiaMedia.SimpleMediumData",
    "category": "type",
    "text": "data = SimpleMediumData(;cp_const=NaN, cv_const=NaN, d_const=NaN, eta_const=NaN,\n                         lambda_const=NaN, a_const=NaN, T0=NaN, MM_const=NaN)\n\nGenerate a SimpleMediumData object containing the data for a SimpleMedium medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasa",
    "page": "Types",
    "title": "ModiaMedia.SingleGasNasa",
    "category": "type",
    "text": "medium = SingleGasNasa(;infos=nothing, constants=nothing, limits=FluidLimits(), data=nothing)\n\nGenerate a SingleGasNasa <: PureSubstance medium object.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.SingleGasNasaData",
    "page": "Types",
    "title": "ModiaMedia.SingleGasNasaData",
    "category": "type",
    "text": "data = SingleGasNasaData(;name=Missing, MM=NaN, Hf=NaN, H0=NaN, Tlimit=NaN,\n                          alow=Missing, blow=Missing, ahigh=Missing,\n                          bhigh=Missing, R=NaN)\n\nGenerate a SingleGasNasaData object containing the data for an ideal Gas based on the NASA Glenn coefficients.\n\n\n\n\n\n"
},

{
    "location": "lib/Types.html#ModiaMedia.Th",
    "page": "Types",
    "title": "ModiaMedia.Th",
    "category": "type",
    "text": "Enumeration Th - defining whether T or h are known as boundary condition\n\n\n\n\n\n"
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
    "text": "state = ThermodynamicState_pT(p,T)\n\nGenerate a ThermodynamicState_pT <: PureSubstanceThermodynamicState object containg pressure p [Pa] and temperature T [K] as states.\n\n\n\n\n\n"
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
    "location": "lib/Functions.html#ModiaMedia.cp_T-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.cp_T",
    "category": "method",
    "text": "cp = cp_T(data::SingleGasNasaData, T) - Compute specific heat capacity at constant pressure from temperature and gas data\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.density",
    "category": "method",
    "text": "d = density(medium,state) - return density for medium from state in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_pTX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_pTX",
    "category": "method",
    "text": "d = density_pTX(medium,p,T,X) - return density for medium from p, T, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_ph-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_ph",
    "category": "method",
    "text": "d = density_ph(medium,p,h) - return density in [kg/m^3] for medium::PureSubstance from pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.density_phX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.density_phX",
    "category": "method",
    "text": "d = density_phX(medium,p,h,X) - return density for medium from p, h, and X or Xi in [kg/m^3]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.h_T",
    "page": "Functions",
    "title": "ModiaMedia.h_T",
    "category": "function",
    "text": "h = h_T(data, T, exclEnthForm=true, refChoice=ReferenceEnthalpy_ZeroAt0K, h_off=0.0)\n\nReturn specific enthalpy from temperature and gas data.\n\nArguments\n\ndata::SingleGasNasaData: Data of the SingleGasNasa medium.\nT::Float64: Temperature in [K].\nexclEnthForm::Bool: If true, enthalpy of formation Hf is not included in specific enthalpy h.\nrefChoice::ReferenceEnthalpy: Choice of reference enthalpy.\nh_off::Float64: User defined offset for reference enthalpy, if referenceEnthalpy = ReferenceEnthalpy_UserDefined\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.pressure",
    "category": "method",
    "text": "p = pressure(medium,state) - return pressure for medium from state in [Pa]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.pressure_dT-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.pressure_dT",
    "category": "method",
    "text": "p = pressure_dT(medium,d,T) - return pressure in [Pa] for medium::PureSubstance from density c in [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.s0_T-Tuple{ModiaMedia.SingleGasNasaData,Float64}",
    "page": "Functions",
    "title": "ModiaMedia.s0_T",
    "category": "method",
    "text": "Compute specific entropy from temperature and gas data\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dT-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dT",
    "category": "method",
    "text": "state = setState_dT(medium, d,T)\n\nGenerate a state object for medium medium::PureSubstance for density d [kg/m^3] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_dTX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_dTX",
    "category": "method",
    "text": "state = setState_dTX(medium, d,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for density d [kg/m^3], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pT-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pT",
    "category": "method",
    "text": "state = setState_pT(medium, p,T)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_pTX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_pTX",
    "category": "method",
    "text": "state = setState_pTX(medium, p,T,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], temperature T [K] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ph-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ph",
    "category": "method",
    "text": "state = setState_ph(medium, p,h)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_phX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_phX",
    "category": "method",
    "text": "state = setState_phX(medium, p,h,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific enthalpy h [J/kg]] and mass fractions vector X or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_ps-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_ps",
    "category": "method",
    "text": "state = setState_ps(medium, p,s)\n\nGenerate a state object for medium medium::PureSubstance for pressure p [Pa] and specific entropy s [J/(kg*K)].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.setState_psX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.setState_psX",
    "category": "method",
    "text": "state = setState_psX(medium, p,s,X)\n\nGenerate a state object for medium medium::AbstractMedium for pressure p [Pa], specific entropy s [J/(kg*K)] and mass fractions vector X  or Xi.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy",
    "category": "method",
    "text": "h = specificEnthalpy(medium,state) - return specific enthalpy for medium from state in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_dT-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy_dT",
    "category": "method",
    "text": "h = specificEnthalpy_dT(medium,d,T) - return specific enthalpy in [J/kg] for medium from density d [kg/m^3] and and temperature T [K].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificEnthalpy_pTX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.specificEnthalpy_pTX",
    "category": "method",
    "text": "h = specificEnthalpy_pTX(medium,p,T,X) - return specific enthalpy for medium from p, T, and X or Xi in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificHeatCapacityCp-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificHeatCapacityCp",
    "category": "method",
    "text": "cp = specificHeatCapacityCp(medium,state) - return specific heat capacity at constant pressure for medium from state in [J/(kg*K)]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.specificInternalEnergy-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.specificInternalEnergy",
    "category": "method",
    "text": "u = specificInternalEnergy(medium,state) - return specific internal energy from medium at state in [J/kg]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardCharacteristics-Tuple{ModiaMedia.AbstractMedium}",
    "page": "Functions",
    "title": "ModiaMedia.standardCharacteristics",
    "category": "method",
    "text": "dict = standardCharacteristics(medium::AbstractMedium)\n\nReturn a dict::dict{AbstractString,Any} dictionary with the most  important characteristics of the medium as vectors.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.standardPlot-Tuple{ModiaMedia.AbstractMedium}",
    "page": "Functions",
    "title": "ModiaMedia.standardPlot",
    "category": "method",
    "text": "standardPlot(medium::AbstractMedium; figure=1)\n\nPlot the standardCharacteristics(medium) of the medium.\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature-Tuple{ModiaMedia.AbstractMedium,ModiaMedia.ThermodynamicState}",
    "page": "Functions",
    "title": "ModiaMedia.temperature",
    "category": "method",
    "text": "T = temperature(medium,state) - return temperature for medium from state in [K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_ph-Tuple{ModiaMedia.PureSubstance,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.temperature_ph",
    "category": "method",
    "text": "T = temperature_ph(medium,p,h) - return temperature in [K] for medium::PureSubstance from pressure p [Pa] and specific enthalpy h [J/kg].\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#ModiaMedia.temperature_phX-Tuple{ModiaMedia.AbstractMedium,Any,Any,Any}",
    "page": "Functions",
    "title": "ModiaMedia.temperature_phX",
    "category": "method",
    "text": "T = temperature_phX(medium,p,h,X) - return temperature for medium from p, h, and X or Xi in [K]\n\n\n\n\n\n"
},

{
    "location": "lib/Functions.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [ModiaMedia]\r\nOrder   = [:function]"
},

]}
