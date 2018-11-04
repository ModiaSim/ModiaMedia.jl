#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

### Data structures that are additionally available for media <: PartialPureSubstance -----------------------

"""
    state = ThermodynamicState_pT(p,T)

Generate a `ThermodynamicState_pT <: PureSubstanceThermodynamicState` object containg
pressure `p` [Pa] and temperature `T` [K] as states.
"""
mutable struct ThermodynamicState_pT <: PureSubstanceThermodynamicState
    p::Float64
    T::Float64
end


""" 
    state = setState_pT(medium, p,T)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and temperature `T` [K].
"""
setState_pT( m::PureSubstance,p,T) = setState_pTX(m,p,T,fill(0.0,0))


""" 
    state = setState_ph(medium, p,h)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and specific enthalpy `h` [J/kg].
"""
setState_ph( m::PureSubstance,p,h) = setState_phX(m,p,h,fill(0.0,0))


""" 
    state = setState_ps(medium, p,s)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and specific entropy `s` [J/(kg*K)].
"""
setState_ps( m::PureSubstance,p,s) = setState_psX(m,p,s,fill(0.0,0))


""" 
    state = setState_dT(medium, d,T)

Generate a state object for medium `medium::PureSubstance` for
density `d` [kg/m^3] and temperature `T` [K].
"""
setState_dT( m::PureSubstance,d,T) = setState_dTX(m,d,T,fill(0.0,0))



### Medium functions that are additionally available for media <: PartialPureSubstance -----------------------

"d = density_ph(medium,p,h) - return density in [kg/m^3] for `medium::PureSubstance` from pressure `p` [Pa] and specific enthalpy `h` [J/kg]."
density_ph(m::PureSubstance, p,h) = density(m, setState_ph(m,p,h))

"T = temperature_ph(medium,p,h) - return temperature in [K] for `medium::PureSubstance` from pressure `p` [Pa] and specific enthalpy `h` [J/kg]."
temperature_ph(m::PureSubstance, p,h) = temperature(m, setState_ph(m,p,h))

"p = pressure_dT(medium,d,T) - return pressure in [Pa] for `medium::PureSubstance` from density `c` in [kg/m^3] and temperature `T` [K]."
pressure_dT(m::PureSubstance, d,T) = pressure(m, setState_dT(m,d,T))

"h = specificEnthalpy_dT(medium,d,T) - return specific enthalpy in [J/kg] for `medium` from density `d` [kg/m^3] and and temperature `T` [K]."
specificEnthalpy_dT(m::PureSubstance, d,T) = specificEnthalpy(m, setState_dT(m,d,T))


