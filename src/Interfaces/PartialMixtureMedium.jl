#
# This file is part of module
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

### Data structures that are additionally available for media <: PartialMixtureMedium -----------------------


"""
    R = gasConstant(state)

Return gas constant of MixtureMedium from
`state::MixtureThermodynamicState` in [J/mol.K]
"""
gasConstant(medium::MixtureMedium, state::MixtureThermodynamicState) = undefinedFunction("gasConstant", state)
gasConstant(                       state::MixtureThermodynamicState) = gasConstant(state.medium,state)

# moleToMassFractions: use inner constructor with fluid information.
# massToMoleFractions: use inner constructor with fluid information from abstract medium, since we have that type.