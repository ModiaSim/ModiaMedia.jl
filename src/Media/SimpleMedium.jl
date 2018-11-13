#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

"""
    data = SimpleMediumData(;cp_const=NaN, cv_const=NaN, d_const=NaN, eta_const=NaN,
                             lambda_const=NaN, a_const=NaN, T0=NaN, MM_const=NaN)

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
    T0::Float64             # "Zero enthalpy temperature";
    MM_const                # "Molar mass";

    SimpleMediumData(;cp_const=NaN, cv_const=NaN, d_const=NaN, eta_const=NaN,
                      lambda_const=NaN, a_const=NaN, T0=NaN, MM_const=NaN) =
         new(cp_const, cv_const, d_const, eta_const, lambda_const, a_const, T0, MM_const)
end



"""
    medium = SimpleMedium(; mediumName     = Missing,
                            reference_p    = 101325,
                            reference_T    = 298.15,
                            p_default      = 101325,
                            T_default      = 293.15,
                            fluidConstants = nothing, 
                            fluidLimits    = FluidLimits(), 
                            data           = nothing)

Generate a `SimpleMedium <: PureSubstance` medium object.
"""
struct SimpleMedium <: PureSubstance
    infos::FluidInfos
    fluidConstants::SVector{1,BasicFluidConstants}
    fluidLimits::FluidLimits
    data::SimpleMediumData

    function SimpleMedium(; mediumName=Missing,
                            reference_p=101325,
                            reference_T=298.15,
                            p_default=101325,
                            T_default=293.15,
                            fluidConstants=nothing, 
                            fluidLimits=FluidLimits(), 
                            data=nothing)

        infos = FluidInfos(mediumName           = mediumName,
                           substanceNames       = [mediumName],
                           extraPropertiesNames = fill("",0),
                           ThermoStates         = IndependentVariables_T,
                           baseProperties       = :BaseProperties_SimpleMedium,
                           singleState          = true,
                           reducedX             = true,
                           fixedX               = false,
                           reference_p          = reference_p,
                           reference_T          = reference_T,
                           reference_X          = fill(1.0,1),
                           p_default            = p_default,
                           T_default            = T_default,
                           h_default            = specificEnthalpy(data, ThermodynamicState_pT(p_default,T_default)))

        fluidLimits.HMIN = specificEnthalpy(data, ThermodynamicState_pT(p_default,fluidLimits.TMIN))
        fluidLimits.HMAX = specificEnthalpy(data, ThermodynamicState_pT(p_default,fluidLimits.TMAX))

        new(infos, fill(fluidConstants,1), fluidLimits, data)
    end
end


setState_pTX(m::SimpleMedium,p,T,X) = ThermodynamicState_pT(p,T)
setState_phX(m::SimpleMedium,p,h,X) = ThermodynamicState_pT(p,m.data.T0+h/m.data.cp_const)
setState_psX(m::SimpleMedium,p,s,X) = ThermodynamicState_pT(p,exp(s/m.data.cp_const + log(m.infos.reference_T)))
setState_dTX(m::SimpleMedium,d,T,X) = error("From setState_dTX: Pressure cannot be computed from temperature and density for the incompressible fluid $(m.infos.mediumName)!")


specificEnthalpy(data::SimpleMediumData, state::ThermodynamicState_pT)::Float64 = data.cp_const*(state.T - data.T0)
pressure(               m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = state.p
temperature(            m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = state.T
density(                m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.d_const
specificEnthalpy(       m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = specificEnthalpy(m.data,state)
specificInternalEnergy( m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.cv_const*(state.T - m.data.T0)
specificHeatCapacityCp( m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const



function standardCharacteristics(m::SimpleMedium)::Dict{AbstractString,Any}
    p = m.infos.reference_p
    T = collect( range(m.fluidLimits.TMIN, stop=m.fluidLimits.TMAX, length=101) )
    h = zeros(length(T))

    for i in 1:length(T)
        h[i] = specificEnthalpy(m, setState_pT(m,p,T[i]))
    end

    mediumDict = Dict{AbstractString,Any}()
    mediumDict["T"] = uconvert.(u"°C", T*1u"K")
    mediumDict["h"] = h*U"J/kg"

    return mediumDict
end


function standardPlot(m::SimpleMedium; figure=1) 
    mediumDict = standardCharacteristics(m)
    ModiaMath.plot(mediumDict, "h", xAxis="T", heading=m.infos.mediumName, figure=figure)
end
