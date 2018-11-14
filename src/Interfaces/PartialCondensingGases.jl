#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

### Data structures that are additionally available for media <: PartialMixtureMedium -----------------------

saturationPressure(m::AbstractMedium, TSat) = undefinedFunction("saturation pressure", m)