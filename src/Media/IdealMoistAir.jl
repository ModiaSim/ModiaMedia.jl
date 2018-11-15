#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#
using LinearAlgebra

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
                           h_default            = specificEnthalpy(SVector{2,SingleGasNasaData}(data), 
                                                                   ThermodynamicState_pTX(p_default, T_default, X_default)))

        fluidLimits = FluidLimits(TMIN = data[1].T_min, 
                                  TMAX = data[1].T_max,
                                  HMIN = specificEnthalpy(SVector{2,SingleGasNasaData}(data), 
                                                          ThermodynamicState_pTX(p_default, data.T_min, X_default)),
                                  HMAX = specificEnthalpy(SVector{2,SingleGasNasaData}(data), 
                                                          ThermodynamicState_pTX(p_default, data.T_max, X_default)))

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

function specificEnthalpy(data::SVector{2,SingleGasNasaData}, state::ThermodynamicState_pTX)
    p = state.p 
    T = state.T 
    X = state.X 

    # steam.MM/dryAir.MM 
    k_mair = data[1].fluidConstants[1].molarMass/data[2].fluidConstants[1].molarMass 

    p_steam_sat = saturationPressure(T)
    X_sat = min(p_steam_sat*k_mair/max(100*eps(), p - p_steam_sat)*(1-X[1]), 1.0)
    X_liquid = max(X[1] - X_sat, 0.0);
    X_steam = X[1] - X_liquid;
    X_air = 1 - X[1];   

    h = dot([h_T(data[1], T, h_off=46479.819 + 2501014.5); 
             h_T(data[2], T, h_off=25104.684)], 
            [X_steam; X_air]) 
        + enthalpyOfWater(T)*X_liquid

    return h 
end

setState_pTX(m::IdealMoistAirMedium,p,T,X) = ThermodynamicState_pTX(p,T,X)
#setState_phX(m::SimpleIdealGasMedium,p,h,X) = ThermodynamicState_pT(p,m.data.T0+h/m.data.cp_const)
#setState_psX(m::SimpleIdealGasMedium,p,s,X) = ThermodynamicState_pT(p,exp(s/m.data.cp_const + log(m.infos.reference_T)+m.infos.R_gas*log(p/m.infos.reference_p)))
#setState_dTX(m::SimpleIdealGasMedium,d,T,X) = ThermodynamicState_pT(d*m.infos.R_gas*T,T)

pressure(               m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = state.p
temperature(            m::IdealMoistAirMedium, state::ThermodynamicState_pTX)::Float64 = state.T
#density(                m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = state.p/(m.data.R_gas*state.T)
#specificEnthalpy(       m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = specificEnthalpy(m.data,state)
#specificInternalEnergy( m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const*(state.T-m.data.T0)-m.data.R_gas*state.T
#specificEntropy(        m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const*log(state.T/m.data.T0)-m.data.R_gas*log(state.p/m.infos.reference_p)
#specificGibbsEnergy(    m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const*(state.T-m.data.T0)-state.T*specificEntropy(m,state)
#specificHelmholtzEnergy(m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = specificInternalEnergy(state)-state.T*specificEntropy(state)
#dynamicViscosity(       m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.eta_const
#thermalConductivity(    m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.lambda_const
#specificHeatCapacityCp( m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const
#specificHeatCapacityCv( m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cv_const
#isentropicExponent(     m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const/m.data.cv_const
#velocityOfSound(        m::SimpleIdealGasMedium, state::ThermodynamicState_pT)::Float64 = sqrt(m.data.cp_const/m.data.cv_const*m.data.R_gas*state.T)


#density_pT(      m::SimpleIdealGasMedium,p,T)::Float64 = p/(m.data.R_gas*T)
#density_pT_der_1(m::SimpleIdealGasMedium,p,T)::Float64 = 0.0
#density_pT_der_2(m::SimpleIdealGasMedium,p,T)::Float64 = 1.0/(m.data.R_gas*T)
#density_pT_der_3(m::SimpleIdealGasMedium,p,T)::Float64 = -p/(m.data.R_gas*T*T)

#specificInternalEnergy_T(      m::SimpleIdealGasMedium,T)::Float64 = m.data.cp_const*(T-m.data.T0)-m.data.R_gas*T
#specificInternalEnergy_T_der_1(m::SimpleIdealGasMedium,T)::Float64 = 0.0
#specificInternalEnergy_T_der_2(m::SimpleIdealGasMedium,T)::Float64 = m.data.cp_const - m.data.R_gas


#function standardCharacteristics(m::SimpleIdealGasMedium)::Dict{AbstractString,Any}
#    p_ref = m.infos.reference_p
#    T     = collect( range(m.fluidLimits.TMIN, stop=min(1600.0, m.fluidLimits.TMAX), length=501) )
#    p     = [0.5e5, 1.0e5, 2.0e5]
#    nT    = length(T)
#    np    = length(p)
#    h     = zeros(nT)
#    u     = zeros(nT)
#    cp    = zeros(nT)
#    cv    = zeros(nT)
#    d     = zeros(nT,np)

#    for i in 1:nT
#        state = setState_pT(m,p_ref,T[i])
#        h[i]  = specificEnthalpy(      m, state)
#        u[i]  = specificInternalEnergy(m, state)
#        cp[i] = specificHeatCapacityCp(m, state)
#        cv[i] = specificHeatCapacityCv(m, state)
#    end

#    for j in 1:np
#        for i in 1:nT
#            d[i,j] = to_DensityDisplayUnit( density(m, ThermodynamicState_pT(p[j],T[i]) ) )     
#        end
#    end       

#    mediumDict = Dict{AbstractString,Any}()
#    mediumDict["T"]  = uconvert.(u"°C", T*1U"K")
#    mediumDict["h"]  = h*1U"J/kg"
#    mediumDict["u"]  = u*1U"J/kg"
#    mediumDict["cp"] = cp*1U"J/(kg*K)"
#    mediumDict["cv"] = cv*1U"J/(kg*K)"
#    mediumDict["d(p=0.5 bar)"] = d[:,1]*1U"g/cm^3"
#    mediumDict["d(p=1.0 bar)"] = d[:,2]*1U"g/cm^3"
#    mediumDict["d(p=2.0 bar)"] = d[:,3]*1U"g/cm^3"
#    return mediumDict
#end


#function standardPlot(m::SimpleIdealGasMedium; figure=1) 
#    mediumDict = standardCharacteristics(m)
#    ModiaMath.plot(mediumDict, [("h", "u"), ("cp","cv"), ("d(p=0.5 bar)" ,
#                                                          "d(p=1.0 bar)" ,
#                                                          "d(p=2.0 bar)")], xAxis="T", heading=m.infos.mediumName, figure=figure)
#end