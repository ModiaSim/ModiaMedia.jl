#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#


"""
    data = SimpleMediumData(;cp_const=nothing, cv_const=nothing, d_const=nothing, eta_const=nothing,
                             lambda_const=nothing, a_const=nothing, T0=nothing, MM_const=nothing)

Generate a `SimpleMediumData` object containing the data
for a SimpleMedium medium.
"""
mutable struct SimpleMediumData
    cp_const::Float64       # "Constant specific heat capacity at constant pressure";
    cv_const::Float64       # "Constant specific heat capacity at constant volume";
    d_const::Float64        # "Constant density";
    eta_const::Float64      # "Constant dynamic viscosity";
    lambda_const::Float64   # "Constant thermal conductivity";
    a_const::Float64        # "Constant velocity of sound";
    T_min::Float64          # Minimum temperature valid for medium model 
    T_max::Float64          # Maximum temperature valid for medium model
    T0::Float64             # "Zero enthalpy temperature";
    MM_const                # "Molar mass";

    SimpleMediumData(;cp_const=nothing, cv_const=nothing, d_const=nothing, eta_const=nothing,
                      lambda_const=nothing, a_const=nothing, T_min=nothing, T_max=nothing, 
                      T0=nothing, MM_const=nothing) =
         new(cp_const, cv_const, d_const, eta_const, lambda_const, a_const, T_min, T_max, T0, MM_const)
end



"""
    medium = SimpleMedium(; mediumName     = nothing,
                            reference_p    = 101325.0,
                            reference_T    = 298.15,
                            p_default      = 101325.0,
                            T_default      = 293.15,
                            fluidConstants = nothing, 
                            data           = nothing)

Generate a `SimpleMedium <: PureSubstance` medium object.
"""
mutable struct SimpleMedium <: PureSubstance
    infos::FluidInfos
    fluidConstants::SVector{1,BasicFluidConstants}
    fluidLimits::FluidLimits
    data::SimpleMediumData
end


mutable struct SimpleMediumState <: ThermodynamicState
    Medium::SimpleMedium
    p::Float64
    T::Float64
end


function SimpleMedium(; mediumName=nothing,
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
                       baseProperties       = :BaseProperties_SimpleMedium,
                       singleState          = true,
                       reducedX             = true,
                       fixedX               = false,
                       reference_p          = reference_p,
                       reference_T          = reference_T,
                       reference_X          = fill(1.0,1),
                       p_default            = p_default,
                       T_default            = T_default)

    fluidLimits = FluidLimits(TMIN = data.T_min, 
                              TMAX = data.T_max)

    Medium = SimpleMedium(infos, fill(fluidConstants,1), fluidLimits, data)

    infos.h_default  = specificEnthalpy(data, SimpleMediumState(Medium,p_default,T_default))
    fluidLimits.HMIN = specificEnthalpy(data, SimpleMediumState(Medium,p_default,data.T_min))
    fluidLimits.HMAX = specificEnthalpy(data, SimpleMediumState(Medium,p_default,data.T_max))

    return Medium
end

setState_pTX(m::SimpleMedium,p,T,X) = SimpleMediumState(m,p,T)
setState_phX(m::SimpleMedium,p,h,X) = SimpleMediumState(m,p,m.data.T0+h/m.data.cp_const)
setState_psX(m::SimpleMedium,p,s,X) = SimpleMediumState(m,p,exp(s/m.data.cp_const + log(m.infos.reference_T)))
setState_dTX(m::SimpleMedium,d,T,X) = error("From setState_dTX: Pressure cannot be computed from temperature and density for the incompressible fluid $(m.infos.mediumName)!")
isenthalpicState(m::SimpleMedium, state::SimpleMediumState, dp::Float64) = SimpleMediumState(m,state.p+dp, state.T)



specificEnthalpy(data::SimpleMediumData, state::SimpleMediumState)::Float64 = data.cp_const*(state.T - data.T0)
pressure(               m::SimpleMedium, state::SimpleMediumState)::Float64 = state.p
temperature(            m::SimpleMedium, state::SimpleMediumState)::Float64 = state.T
density(                m::SimpleMedium, state::SimpleMediumState)::Float64 = m.data.d_const
specificEnthalpy(       m::SimpleMedium, state::SimpleMediumState)::Float64 = specificEnthalpy(m.data,state)
specificInternalEnergy( m::SimpleMedium, state::SimpleMediumState)::Float64 = m.data.cv_const*(state.T - m.data.T0)
specificHeatCapacityCp( m::SimpleMedium, state::SimpleMediumState)::Float64 = m.data.cp_const
dynamicViscosity(       m::SimpleMedium, state::SimpleMediumState)::Float64 = m.data.eta_const

specificEnthalpy_T(m::SimpleMedium, T)::Float64 = m.data.cp_const*(T - m.data.T0)


### Functions for mass and energy balance

density(      m::SimpleMedium)::Float64 = m.data.d_const
density_der_1(m::SimpleMedium)::Float64 = 0.0

specificInternalEnergy_T(      m::SimpleMedium,T)::Float64 = m.data.cv_const*(T - m.data.T0)
specificInternalEnergy_T_der_1(m::SimpleMedium,T)::Float64 = 0.0
specificInternalEnergy_T_der_2(m::SimpleMedium,T)::Float64 = m.data.cv_const


function standardCharacteristics(m::SimpleMedium)::Dict{AbstractString,Any}
    p = m.infos.reference_p
    T = collect( range(m.fluidLimits.TMIN, stop=m.fluidLimits.TMAX, length=101) )
    h = zeros(length(T))

    for i in 1:length(T)
        h[i] = specificEnthalpy(setState_pT(m,p,T[i]))
    end

    mediumDict = Dict{AbstractString,Any}()
    mediumDict["T"] = uconvert.(u"°C", T*1u"K")
    mediumDict["h"] = h*u"J/kg"

    return mediumDict
end


function standardPlot(m::SimpleMedium; figure=1) 
    mediumDict = standardCharacteristics(m)
    ModiaMath.plot(mediumDict, "h", xAxis="T", heading=m.infos.mediumName, figure=figure)
end
