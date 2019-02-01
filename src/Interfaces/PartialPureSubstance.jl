#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

### Data structures that are additionally available for media <: PartialPureSubstance -----------------------

const X_dummy = fill(0.0,0)


""" 
    state = setState_pT(medium, p,T)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and temperature `T` [K].
"""
setState_pT( m::PureSubstance,p,T) = setState_pTX(m,p,T,X_dummy)



""" 
    setState_pT!(state, p,T)

Update the `state::ThermodynamicState` of a `PureSubstance` medium with 
pressure `p` [Pa] and temperature `T` [K].
"""
setState_pT!(state::ThermodynamicState,p,T) = setState_pTX!(state,p,T,X_dummy)



""" 
    state = setState_ph(medium, p,h)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and specific enthalpy `h` [J/kg].
"""
setState_ph(m::PureSubstance,p,h) = setState_phX(m,p,h,X_dummy)


""" 
    setState_ph!(state, p,h)

Update the `state::ThermodynamicState` of a `PureSubstance` medium with 
pressure `p` [Pa] and specific enthalpy `h` [J/kg]].
"""
setState_ph!(state::ThermodynamicState,p,h) = setState_phX!(state,p,h,X_dummy)



""" 
    state = setState_ps(medium, p,s)

Generate a state object for medium `medium::PureSubstance` for
pressure `p` [Pa] and specific entropy `s` [J/(kg*K)].
"""
setState_ps( m::PureSubstance,p,s) = setState_psX(m,p,s,X_dummy)



""" 
    setState_ps!(state, p,s)

Update the `state::ThermodynamicState` of a `PureSubstance` medium with 
pressure `p` [Pa] and specific entropy `s` [J/(kg*K)].
"""
setState_ps!(state::ThermodynamicState,p,s) = setState_psX!(state,p,s,X_dummy)



""" 
    state = setState_dT(medium, d,T)

Generate a state object for medium `medium::PureSubstance` for
density `d` [kg/m^3] and temperature `T` [K].
"""
setState_dT( m::PureSubstance,d,T) = setState_dTX(m,d,T,X_dummy)


""" 
    setState_dT!(state, d,T)

Update the `state::ThermodynamicState` of a `PureSubstance` medium with 
density `d` [kg/m^3] and temperature `T` [K].
"""
setState_dT!(state::ThermodynamicState,d,T) = setState_dTX!(state,d,T,X_dummy)




### Medium functions that are additionally available for media <: PartialPureSubstance -----------------------

"d = density_ph(medium,p,h) - return density in [kg/m^3] for `medium::PureSubstance` from pressure `p` [Pa] and specific enthalpy `h` [J/kg]."
density_ph(m::PureSubstance, p,h) = density(setState_ph(m,p,h))

"T = temperature_ph(medium,p,h) - return temperature in [K] for `medium::PureSubstance` from pressure `p` [Pa] and specific enthalpy `h` [J/kg]."
temperature_ph(m::PureSubstance, p,h) = temperature(setState_ph(m,p,h))

"p = pressure_dT(medium,d,T) - return pressure in [Pa] for `medium::PureSubstance` from density `c` in [kg/m^3] and temperature `T` [K]."
pressure_dT(m::PureSubstance, d,T) = pressure(setState_dT(m,d,T))

"h = specificEnthalpy_dT(medium,d,T) - return specific enthalpy in [J/kg] for `medium` from density `d` [kg/m^3] and and temperature `T` [K]."
specificEnthalpy_dT(m::PureSubstance, d,T) = specificEnthalpy(setState_dT(m,d,T))




### Derivatives needed for mass and energy balance ------------------------------------------------------

density_ph_der_1(m::PureSubstance, p, h) = 0.0
density_ph_der_2(m::PureSubstance, p, h) = undefinedFunction("density_ph_der_2", m)
density_ph_der_3(m::PureSubstance, p, h) = undefinedFunction("density_ph_der_3", m)



