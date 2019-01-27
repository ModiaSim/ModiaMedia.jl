#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

### Data structures that are additionally available for media <: PartialMixtureMedium -----------------------

"""
    state = ThermodynamicState_pTX(p,T,X)

Generate a `ThermodynamicState_pT <: MixtureThermodynamicState` object containing
pressure `p` [Pa], temperature `T` [K], and a vector of mass fractions `X` as states.
"""
mutable struct ThermodynamicState_pTX <: MixtureThermodynamicState
    Medium::AbstractMedium
    p::Float64
    T::Float64
    X::Vector{Float64}
end

"R = gas_constant(medium,state) - return gas constant for `medium` from `state` in [J/mol.K]"
# should this be MixtureMedium?
# should this take MixtureThermodynamicState? 
gas_constant(m::AbstractMedium, state::ThermodynamicState) = undefinedFunction("gas constant", m)
gas_constant(                   state::ThermodynamicState) = undefinedFunction("gas constant", state)

# moleToMassFractions: use inner constructor with fluid information. 
# massToMoleFractions: use inner constructor with fluid information from abstract medium, since we have that type.