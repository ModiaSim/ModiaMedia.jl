#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

"""
    data = SimpleIdealGasMediumData(;cp_const=nothing, R_gas=nothing, MM_const=nothing, 
                                     eta_const=nothing, lambda_const=nothing, T_min=nothing,
                                     T_max=nothing, T0=nothing)

Generate a `SimpleIdealGasMediumData` object containing the data
for a SimpleIdealGas medium.
"""
mutable struct SimpleIdealGasMediumData
    cp_const::Float64       # Constant specific heat capacity at constant pressure
    cv_const::Float64       # Constant specific heat capacity at constant volume
    R_gas::Float64          # Medium specific gas constant
    MM_const::Float64       # Molar mass
    eta_const::Float64      # Constant dynamic viscosity
    lambda_const::Float64   # Constant thermal conductivity
    T_min::Float64          # Minimum temperature valid for medium model 
    T_max::Float64          # Maximum temperature valid for medium model
    T0::Float64             # Zero enthalpy temperature

    function SimpleIdealGasMediumData(;cp_const=nothing, R_gas=nothing, MM_const=nothing, 
                                       eta_const=nothing, lambda_const=nothing, T_min=nothing, T_max=nothing, T0=nothing)
         cv_const = cp_const - R_gas
         new(cp_const, cv_const, R_gas, MM_const, eta_const, lambda_const, 
             T_min, T_max, T0)
    end
end



"""
    medium = SimpleIdealGasMedium(; mediumName     = nothing,
                                    reference_p    = 101325.0,
                                    reference_T    = 298.15,
                                    p_default      = 101325.0,
                                    T_default      = 293.15,
                                    fluidConstants = nothing, 
                                    data           = nothing)

Generate a `SimpleIdealGasMedium <: PureSubstance` medium object.
"""
mutable struct SimpleIdealGasMedium <: PureSubstance
    infos::FluidInfos
    fluidConstants::SVector{1,BasicFluidConstants}
    fluidLimits::FluidLimits
    data::SimpleIdealGasMediumData
end


"""
    state = SimpleIdealGasMediumState(Medium, p, T)

Generate a `SimpleIdealGasMediumState <: ThermodynamicState` object containing
pressure `p` [Pa] and temperature `T` [K] as thermodynamic states.
"""
mutable struct SimpleIdealGasMediumState <: ThermodynamicState
    Medium::SimpleIdealGasMedium
    p::Float64
    T::Float64
end


function SimpleIdealGasMedium(; mediumName=nothing,
                                reference_p=101325.0,
                                reference_T=298.15,
                                p_default=101325.0,
                                T_default=293.15,
                                fluidConstants=nothing, 
                                data=nothing)

    infos = FluidInfos(mediumName           = mediumName,
                       substanceNames       = [mediumName],
                       extraPropertiesNames = fill("",0),
                       ThermoStates         = IndependentVariables_pT,
                       singleState          = false,
                       reducedX             = true,
                       fixedX               = false,
                       reference_p          = reference_p,
                       reference_T          = reference_T,
                       reference_X          = fill(1.0,1),
                       p_default            = p_default,
                       T_default            = T_default)

    fluidLimits = FluidLimits(TMIN = data.T_min, 
                              TMAX = data.T_max)

    Medium = SimpleIdealGasMedium(infos, fill(fluidConstants,1), fluidLimits, data)


    infos.h_default  = specificEnthalpy(data, SimpleIdealGasMediumState(Medium,p_default,T_default))
    fluidLimits.HMIN = specificEnthalpy(data, SimpleIdealGasMediumState(Medium,p_default,data.T_min))
    fluidLimits.HMAX = specificEnthalpy(data, SimpleIdealGasMediumState(Medium,p_default,data.T_max))

    return Medium
end



setState_pTX(m::SimpleIdealGasMedium,p,T,X) = SimpleIdealGasMediumState(m,p,T)
setState_phX(m::SimpleIdealGasMedium,p,h,X) = SimpleIdealGasMediumState(m,p,m.data.T0+h/m.data.cp_const)
setState_psX(m::SimpleIdealGasMedium,p,s,X) = SimpleIdealGasMediumState(m,p,exp(s/m.data.cp_const + log(m.infos.reference_T)+m.infos.R_gas*log(p/m.infos.reference_p)))
setState_dTX(m::SimpleIdealGasMedium,d,T,X) = SimpleIdealGasMediumState(m,d*m.data.R_gas*T,T)
isenthalpicState(m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState, dp::Float64) = SimpleIdealGasMediumState(m, state.p+dp, state.T)


setState_pTX!(state::SimpleIdealGasMediumState,p,T,X) = begin state.p=p; state.T=T; nothing end
setState_phX!(state::SimpleIdealGasMediumState,p,h,X) = begin state.p=p; state.T=state.Medium.data.T0+h/state.Medium.data.cp_const; nothing end
setState_psX!(state::SimpleIdealGasMediumState,p,s,X) = begin state.p=p; state.T=exp(s/state.Medium.data.cp_const + log(state.Medium.infos.reference_T)); nothing end
setState_dTX!(state::SimpleIdealGasMediumState,d,T,X) = begin state.p=d*state.Medium.data.R_gas*T; state.T=T; nothing end
isenthalpicState!(state_b::SimpleIdealGasMediumState, state_a::SimpleIdealGasMediumState, dp::Float64) = begin state_b.p = state_a.p+dp; state_b.T = state_a.T; nothing end


specificEnthalpy(data::SimpleIdealGasMediumData, state::SimpleIdealGasMediumState)::Float64 = data.cp_const*(state.T - data.T0)
pressure(               m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = state.p
temperature(            m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = state.T
density(                m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = state.p/(m.data.R_gas*state.T)
specificEnthalpy(       m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = specificEnthalpy(m.data,state)
specificInternalEnergy( m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cp_const*(state.T-m.data.T0)-m.data.R_gas*state.T
specificEntropy(        m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cp_const*log(state.T/m.data.T0)-m.data.R_gas*log(state.p/m.infos.reference_p)
specificGibbsEnergy(    m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cp_const*(state.T-m.data.T0)-state.T*specificEntropy(m,state)
specificHelmholtzEnergy(m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = specificInternalEnergy(state)-state.T*specificEntropy(state)
dynamicViscosity(       m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.eta_const
thermalConductivity(    m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.lambda_const
specificHeatCapacityCp( m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cp_const
specificHeatCapacityCv( m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cv_const
isentropicExponent(     m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = m.data.cp_const/m.data.cv_const
velocityOfSound(        m::SimpleIdealGasMedium, state::SimpleIdealGasMediumState)::Float64 = sqrt(m.data.cp_const/m.data.cv_const*m.data.R_gas*state.T)


density_pT(      m::SimpleIdealGasMedium,p,T)::Float64 = p/(m.data.R_gas*T)
density_pT_der_1(m::SimpleIdealGasMedium,p,T)::Float64 = 0.0
density_pT_der_2(m::SimpleIdealGasMedium,p,T)::Float64 = 1.0/(m.data.R_gas*T)
density_pT_der_3(m::SimpleIdealGasMedium,p,T)::Float64 = -p/(m.data.R_gas*T*T)

specificInternalEnergy_T(      m::SimpleIdealGasMedium,T)::Float64 = m.data.cp_const*(T-m.data.T0)-m.data.R_gas*T
specificInternalEnergy_T_der_1(m::SimpleIdealGasMedium,T)::Float64 = 0.0
specificInternalEnergy_T_der_2(m::SimpleIdealGasMedium,T)::Float64 = m.data.cp_const - m.data.R_gas

specificEnthalpy_T(m::SimpleIdealGasMedium, T)::Float64 = m.data.cp_const*(T - m.data.T0)


function standardCharacteristics(m::SimpleIdealGasMedium)::Dict{AbstractString,Any}
    p_ref = m.infos.reference_p
    T     = collect( range(m.fluidLimits.TMIN, stop=min(1600.0, m.fluidLimits.TMAX), length=501) )
    p     = [0.5e5, 1.0e5, 2.0e5]
    nT    = length(T)
    np    = length(p)
    h     = zeros(nT)
    u     = zeros(nT)
    cp    = zeros(nT)
    cv    = zeros(nT)
    d     = zeros(nT,np)

    for i in 1:nT
        state = setState_pT(m,p_ref,T[i])
        h[i]  = specificEnthalpy(      m, state)
        u[i]  = specificInternalEnergy(m, state)
        cp[i] = specificHeatCapacityCp(m, state)
        cv[i] = specificHeatCapacityCv(m, state)
    end

    for j in 1:np
        for i in 1:nT
            d[i,j] = to_DensityDisplayUnit( density(SimpleIdealGasMediumState(m,p[j],T[i]) ) )     
        end
    end       

    mediumDict = Dict{AbstractString,Any}()
    mediumDict["T"]  = uconvert.(u"°C", T*1u"K")
    mediumDict["h"]  = h*1u"J/kg"
    mediumDict["u"]  = u*1u"J/kg"
    mediumDict["cp"] = cp*1u"J/(kg*K)"
    mediumDict["cv"] = cv*1u"J/(kg*K)"
    mediumDict["d(p=0.5 bar)"] = d[:,1]*1u"g/cm^3"
    mediumDict["d(p=1.0 bar)"] = d[:,2]*1u"g/cm^3"
    mediumDict["d(p=2.0 bar)"] = d[:,3]*1u"g/cm^3"
    return mediumDict
end


function standardPlot(m::SimpleIdealGasMedium; figure=1) 
    mediumDict = standardCharacteristics(m)
    ModiaMath.plot(mediumDict, [("h", "u"), ("cp","cv"), ("d(p=0.5 bar)" ,
                                                          "d(p=1.0 bar)" ,
                                                          "d(p=2.0 bar)")], xAxis="T", heading=m.infos.mediumName, figure=figure)
end