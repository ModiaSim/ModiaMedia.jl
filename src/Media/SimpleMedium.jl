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
    medium = SimpleMedium(;infos=nothing, constants=nothing, limits=FluidLimits(), data=nothing)

Generate a `SimpleMedium <: PureSubstance` medium object.
"""
struct SimpleMedium <: PureSubstance
    infos::FluidInfos
    constants::SVector{1,BasicFluidConstants}
    limits::FluidLimits
    data::SimpleMediumData

    SimpleMedium(;infos=nothing, constants=nothing, limits=FluidLimits(), data=nothing) =
        new(infos, fill(constants,1), limits, data)
end


setState_pTX(m::SimpleMedium,p,T,X) = ThermodynamicState_pT(p,T)
setState_phX(m::SimpleMedium,p,h,X) = ThermodynamicState_pT(p,m.data.T0+h/m.data.cp_const)
setState_psX(m::SimpleMedium,p,s,X) = ThermodynamicState_pT(p,exp(s/m.data.cp_const + log(m.infos.reference_T)))
setState_dTX(m::SimpleMedium,d,T,X) = error("From setState_dTX: Pressure cannot be computed from temperature and density for the incompressible fluid $(m.infos.mediumName)!")

pressure(              m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = state.p
temperature(           m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = state.T
density(               m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.d_const
specificEnthalpy(      m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const*(state.T - m.data.T0)
specificInternalEnergy(m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.cv_const*(state.T - m.data.T0)
specificHeatCapacityCp(m::SimpleMedium, state::ThermodynamicState_pT)::Float64 = m.data.cp_const


function standardCharacteristics(m::SimpleMedium)::Dict{AbstractString,Any}
    p = m.infos.reference_p
    T = collect( range(m.limits.TMIN, stop=m.limits.TMAX, length=101) )
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
