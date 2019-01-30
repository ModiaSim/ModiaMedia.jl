#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#
using LinearAlgebra

const IdealMoistAir_excludeEnthalpyOfFormation = true                            # If true, enthalpy of formation Hf is not included in specific enthalpy hc
const IdealMoistAir_referenceEnthalpy          = ReferenceEnthalpy_UserDefined   # Choice of reference enthalpy
const IdealMoistAir_h_offset                   = 0.0                             # User defined offset for reference enthalpy, if ReferenceEnthalpy = ReferenceEnthalpy_UserDefined




### Data structures

"""
    data = IdealMoistAirMediumData()

Generate an `IdealMoistAirMediumData` object containing the data
for an moist air model.
"""
mutable struct IdealMoistAirMediumData
    Water::Int                  # "Index of water (in substanceNames, massFractions X, etc.)"
    Air::Int                    # "Index of air (in substanceNames, massFractions X, etc.)";
    k_mair::Float64             # "Ratio of molar weights";
    dryair::SingleGasNasaData   # Data of dry air             
    steam::SingleGasNasaData    # Data of steam
    MMX::SVector{2,Float64}     # "Molar masses of components";

    function IdealMoistAirMediumData()
        dryair = mediumDict["Air"].data
        steam  = mediumDict["H2O"].data

        new(1,2,steam.MM/dryair.MM,dryair, steam, SVector{2,Float64}(steam.MM,dryair.MM))
    end 
end


"""
    medium = IdealMoistAirMedium()

Generate an `IdealMoistAirMedium <: CondensingGases` medium object.
Valid for T = 190 K ... 647 K. 
"""
struct IdealMoistAirMedium <: CondensingGases
    infos::FluidInfos
    fluidConstants::SVector{2,IdealGasFluidConstants}
    fluidLimits::FluidLimits
    data::IdealMoistAirMediumData
end



"""
    state = IdealMoistAirMediumState(Medium, p, T, X)

Generate an `IdealMoistAirMediumState <: MixtureThermodynamicState` object containing
pressure `p` [Pa], temperature `T` [K], and a vector of mass fractions `X[2]` as states.
X[1] is the mass fraction of Steam (Water) and X[2] is the mass fraction of Air.
"""
mutable struct IdealMoistAirMediumState <: MixtureThermodynamicState
    Medium::IdealMoistAirMedium
    p::Float64
    T::Float64
    X::SVector{2,Float64}

    IdealMoistAirMediumState(Medium::IdealMoistAirMedium, p, T, X) = 
        new(Medium, p, T, SVector{2,Float64}(X))
end


function IdealMoistAirMedium()
    infos = FluidInfos(mediumName           = "MoistAir",
                       substanceNames       = ["water", "air"],
                       extraPropertiesNames = fill("",0),
                       ThermoStates         = IndependentVariables_pTX,
                       singleState          = false,
                       reducedX             = true,
                       fixedX               = false,
                       reference_p          = 101325.0,
                       reference_T          = 298.15,
                       reference_X          = [0.01, 0.99],
                       p_default            = 101325.0,
                       T_default            = 293.15,
                       X_default            = [0.01, 0.99])

    # these limits should be fixed after the specific enthalpy can be calculated below 273.15K; the high end also has problems as well. 
    fluidLimits = FluidLimits(TMIN = 190, 
                              TMAX = 647)
    infos.h_default = h_pTX(data, infos.p_default, infos.T_default, infos.X_default)

    # HMIN = h_pTX(IdealMoistAirMediumData(data),p_default,280.0,X_default),
    # HMAX = h_pTX(IdealMoistAirMediumData(data),p_default,600.0,X_default))

    Medium = IdealMoistAirMedium(infos, 
                                 SVector{2,IdealGasFluidConstants}([mediumDict["H2O"].fluidConstants[1], 
                                                                    mediumDict["Air"].fluidConstants[1]]), 
                                 fluidLimits, 
                                 IdealMoistAirMediumData())

    return Medium
end



"Spline interpolation of two functions"
function spliceFunction(
    pos::Float64,                    # Returned value for x-deltax >= 0
    neg::Float64,                    # Returned value for x+deltax <= 0
    x::Float64,                      # Function argument
    deltax::Float64                  # Region around x with spline interpolation
    )::Float64

    scaledX1 = x/deltax
    scaledX = scaledX1*asin(1)
    if scaledX1<=-0.999999999
        y = 0.0
    elseif scaledX1>=0.999999999
        y = 1.0
    else
        y = (tanh(tan(scaledX))+1.0)/2
    end
    return pos*y+(1.0-y)*neg
end


"Return saturation pressure of water as a function of temperature T in the range of 273.16 to 647.096 K"
function saturationPressureLiquid(Tsat::Float64)::Float64
    TCritical=647.096
    pCritical=22.064e6
    r1=(1 - TSat/TCritical)

    a = SVector{6,Float64}(-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502)
    n = SVector{6,Float64}(1.0,1.5,3.0,3.5,4.0,7.5)
    pSat = exp(((a[1]*r1^n[1] + a[2]*r1^n[2] + a[3]*r1^n[3] + a[4]*r1^n[4]
     + a[5]*r1^n[5] + a[6]*r1^n[6])*TCritical)/TSat)*pCritical
    
     return pSat
end


"Return sublimation pressure of water as a function of temperature T between 190 and 273.16 K"
function sublimationPressureIce(Tsat::Float64)::Float64
    Ttriple=273.16                # Triple point temperature
    ptriple=611.657                # Triple point pressure
    r1=Tsat/Ttriple                     # Common subexpression
    a1=-13.9281690
    a2=34.7078238
    n1=-1.5
    n2=-1.25
    psat = exp(a1-a1*r1^n1+a2-a2*r1^n2)*ptriple
    
    return psat
end


 "Return saturation pressure of water as a function of temperature T between 190 and 647.096 K"
saturationPressure(TSat::Float64) = spliceFunction(saturationPressureLiquid(Tsat), sublimationPressureIce(Tsat), Tsat-273.16, 1.0)



"Return specific enthalpy of water (solid/liquid) near atmospheric pressure from temperature T"
enthalpyOfWater(T::Float64) = spliceFunction(4200.0*(T-273.15), 2050.0*(T-273.15)-333000.0, T-273.16, 0.1)



function h_pTX(data::IdealMoistAirMediumData,p::Float64,T::Float64,X::AbstractVector)
    # steam.MM/dryAir.MM 
    k_mair = data.k_mair 

    p_steam_sat = saturationPressure(T)
    X_sat = min(p_steam_sat*k_mair/max(100*eps(), p-p_steam_sat)*(1-X[1]), 1.0)
    X_liquid = max(X[1]-X_sat, 0.0);
    X_steam = X[1]-X_liquid;
    X_air = 1-X[1];   

    h = h_Tlow(data.steam , T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)*X_steam +
        h_Tlow(data.dryair, T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)*X_air +
        enthalpyOfWater(T)*X_liquid

    return h 
end



function u_pTX(data::IdealMoistAirMediumData,p::Float64,T::Float64,X::AbstractVector)
    # steam.MM/dryAir.MM 
    k_mair = data.k_mair

    p_steam_sat = saturationPressure(T)
    X_sat = min(p_steam_sat*k_mair/max(100*eps(), p-p_steam_sat)*(1-X[1]), 1.0)
    X_liquid = max(X[1]-X_sat, 0.0);
    X_steam = X[1]-X_liquid;
    X_air = 1-X[1];       
    R_gas = data.dryair.R*X_air/(1-X_liquid) + data.steam.R*X_steam/(1-X_liquid);

    u = h_Tlow(data.steam , T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)*X_steam +
        h_Tlow(data.dryair, T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)*X_air +
        enthalpyOfWater(T)*X_liquid - R_gas*T
    return u
end


T_phX(data::IdealMoistAirMediumData, p::Float64, h::Float64, X) = ModiaMath.solveOneNonlinearEquation(T->h-h_pTX(data,p,T,X), 190.0, 647; u_nominal=300.0)
p_dTX(data::IdealMoistAirMediumData,d,T,X) = d*(data.steam.R*X[1]+data.dryair.R*(1.0-X[1]))*T

setState_pTX(m::IdealMoistAirMedium,p,T,X) = IdealMoistAirMediumState(m,p,T,X)
setState_phX(m::IdealMoistAirMedium,p,h,X) = IdealMoistAirMediumState(m,p,T_phX(m.data,p,h,X),X)
setState_dTX(m::IdealMoistAirMedium,d,T,X) = IdealMoistAirMediumState(m,p_dTX(m.data,d,T,X),T,X)

pressure(               m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = state.p
temperature(            m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = state.T
gasConstant(            m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = m.data.dryair.R*(1-state.X[1]) + m.data.steam.R*state.X[1]
density(                m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = state.p/(gasConstant(m,state)*state.T)
specificEnthalpy(       m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = h_pTX(m.data,state.p,state.T,state.X)
specificInternalEnergy( m::IdealMoistAirMedium, state::IdealMoistAirMediumState)::Float64 = u_pTX(m.data,state.p,state.T,state.X)

#=
            function dynamicViscosity(fluid::Any, 
                )
                    eta = 1e-6*Polynomials_Temp.evaluateWithRange([9.7391102886305869E-15,-3.1353724870333906E-11,4.3004876595642225E-08,-3.8228016291758240E-05,5.0427874367180762E-02,1.7239260139242528E+01], Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T))
                
                return 
                #= @extends dynamicViscosity
                 =#
            end

=#