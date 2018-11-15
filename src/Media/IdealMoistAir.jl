#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#
using LinearAlgebra

const excludeEnthalpyOfFormation = true                         # If true, enthalpy of formation Hf is not included in specific enthalpy hc
const referenceEnthalpy          = ReferenceEnthalpy_UserDefined   # Choice of reference enthalpy
const h_offset                   = 0.0                          # User defined offset for reference enthalpy, if ReferenceEnthalpy = ReferenceEnthalpy_UserDefined


"""
    medium = SimpleIdealGasMedium(; mediumName     = nothing,
                                    reference_p    = 101325,
                                    reference_T    = 298.15,
                                    p_default      = 101325,
                                    T_default      = 293.15,
                                    fluidConstants = nothing, 
                                    data           = nothing)

Generate a `SimpleIdealGasMedium <: PureSubstance` medium object.
"""

# what is the deal with the baseProperties symbol? 
struct IdealMoistAirMedium <: CondensingGases
    infos::FluidInfos
    fluidConstants::SVector{2,IdealGasFluidConstants}
    fluidLimits::FluidLimits
    data::SVector{2,SingleGasNasaData}

    function IdealMoistAirMedium(; mediumName="MoistAir",
                                    substanceNames=["water", "air"],
                                    reference_p=101325,
                                    reference_T=298.15,
                                    reference_X=[0.01, 0.99],
                                    p_default=101325,
                                    T_default=293.15,
                                    X_default=[0.01, 0.99],
                                    fluidConstants=[mediumDict["H2O"].fluidConstants[1], 
                                                    mediumDict["Air"].fluidConstants[1]], 
                                    data=[mediumDict["H2O"].data, 
                                          mediumDict["Air"].data])

        infos = FluidInfos(mediumName           = mediumName,
                           substanceNames       = substanceNames,
                           extraPropertiesNames = fill("",0),
                           ThermoStates         = IndependentVariables_pTX,
                           baseProperties       = :BaseProperties_SingleGasNasa,
                           singleState          = false,
                           reducedX             = true,
                           fixedX               = false,
                           reference_p          = reference_p,
                           reference_T          = reference_T,
                           reference_X          = reference_X,
                           p_default            = p_default,
                           T_default            = T_default,
                           X_default            = X_default,
                           h_default            = h_pTX(SVector{2,SingleGasNasaData}(data),p_default,T_default,X_default))

        # these limits should be fixed after the specific enthalpy can be calculated below 273.15K; the high end also has problems as well. 
        fluidLimits = FluidLimits(TMIN = max(mediumDict["H2O"].fluidLimits.TMIN, mediumDict["Air"].fluidLimits.TMIN), 
                                  TMAX = min(mediumDict["H2O"].fluidLimits.TMAX, mediumDict["Air"].fluidLimits.TMAX),
                                  HMIN = h_pTX(SVector{2,SingleGasNasaData}(data),p_default,280.0,X_default),
                                  HMAX = h_pTX(SVector{2,SingleGasNasaData}(data),p_default,600.0,X_default))

        new(infos, SVector{2,IdealGasFluidConstants}(fluidConstants), fluidLimits, SVector{2,SingleGasNasaData}(data))
    end
end


function saturationPressure(TSat)
    # this is only valid to 273.15 K (need to add the ice saturation pressure and spline both together)
    TCritical=647.096
    pCritical=22.064e6
    r1=(1 - TSat/TCritical)

    a = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502] 
    n = [1.0,1.5,3.0,3.5,4.0,7.5]
    pSat = exp(((a[1]*r1^n[1] + a[2]*r1^n[2] + a[3]*r1^n[3] + a[4]*r1^n[4]
     + a[5]*r1^n[5] + a[6]*r1^n[6])*TCritical)/TSat)*pCritical
    
     return pSat
end

# enthalpy of water only for liquid condition; need to splice with the enthalpy of ice to get number valid over whole range
enthalpyOfWater(T) = 4200.0

function h_pTX(data::SVector{2,SingleGasNasaData},p,T,X)
    # steam.MM/dryAir.MM 
    k_mair = data[1].MM/data[2].MM 

    p_steam_sat = saturationPressure(T)
    X_sat = min(p_steam_sat*k_mair/max(100*eps(), p-p_steam_sat)*(1-X[1]), 1.0)
    X_liquid = max(X[1]-X_sat, 0.0);
    X_steam = X[1]-X_liquid;
    X_air = 1-X[1];   
    
    println("X_sat=$X_sat\n")
    println("X_liquid=$X_liquid\n")
    println("X_steam=$X_steam\n")
    println("X_air=$X_air\n")


    h = dot([h_T(data[1], T, excludeEnthalpyOfFormation, referenceEnthalpy, (46479.819 + 2501014.5)); 
             h_T(data[2], T, excludeEnthalpyOfFormation, referenceEnthalpy, 25104.684)], 
            [X_steam; X_air]) 
        + enthalpyOfWater(T)*X_liquid

    return h 
end

function u_pTX(data::SVector{2,SingleGasNasaData},p,T,X)
    # steam.MM/dryAir.MM 
    k_mair = data[1].MM/data[2].MM 

    p_steam_sat = saturationPressure(T)
    X_sat = min(p_steam_sat*k_mair/max(100*eps(), p-p_steam_sat)*(1-X[1]), 1.0)
    X_liquid = max(X[1]-X_sat, 0.0);
    X_steam = X[1]-X_liquid;
    X_air = 1-X[1];       
    R_gas = data[2].R*X_air/(1-X_liquid) + data[1].R*X_steam/(1-X_liquid);

    u = h_pTX(data,p,T,X) - R_gas*T
    return u
end


p_dTX(data::SVector{2,SingleGasNasaData},d,T,X) = d*dot([data[1].R; data[2].R], X)*T

setState_pTX(m::IdealMoistAirMedium,p,T,X) = ThermodynamicState_pTX(p,T,X)
setState_dTX(m::IdealMoistAirMedium,d,T,X) = ThermodynamicState_pTX(p_dTX(m.data,d,T,X),T,X)

pressure(               m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = state.p
temperature(            m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = state.T
gasConstant(            m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = m.data[2].R*(1-state.X[1]) + m.data[1].R*state.X[1]
density(                m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = state.p/(gasConstant(m,state)*state.T)
specificEnthalpy(       m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = h_pTX(m.data,state.p,state.T,state.X)
specificInternalEnergy( m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = u_pTX(m.data,state.p,state.T,state.X)
