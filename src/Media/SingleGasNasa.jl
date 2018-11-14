#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

# Dictionary of medium data


const excludeEnthalpyOfFormation = true                         # If true, enthalpy of formation Hf is not included in specific enthalpy hc
const referenceEnthalpy          = ReferenceEnthalpy_ZeroAt0K   # Choice of reference enthalpy
const h_offset                   = 0.0                          # User defined offset for reference enthalpy, if ReferenceEnthalpy = ReferenceEnthalpy_UserDefined




### Data structures

"""
    data = SingleGasNasaData(;name=Missing, MM=nothing, Hf=nothing, H0=nothing, Tlimit=nothing,
                              alow=Missing, blow=Missing, ahigh=Missing,
                              bhigh=Missing, R=nothing)

Generate a `SingleGasNasaData` object containing the data
for an ideal Gas based on the NASA Glenn coefficients.
"""
mutable struct SingleGasNasaData
    name::AbstractString       # Name of ideal gas
    MM::Float64                # SI.MolarMass "Molar mass"
    Hf::Float64                # SI.SpecificEnthalpy "Enthalpy of formation at 298.15K"
    H0::Float64                # SI.SpecificEnthalpy "H0(298.15K) - H0(0K)"
    Tlimit::Float64            # SI.Temperature Tlimit "Temperature limit between low and high data sets"
    alow::SVector{7,Float64}   # "Low temperature coefficients a"
    blow::SVector{2,Float64}   # "Low temperature constants b"
    ahigh::SVector{7,Float64}  # "High temperature coefficients a"
    bhigh::SVector{2,Float64}  # "High temperature constants b"
    R::Float64                 # SI.SpecificHeatCapacity R "Gas constant"

    SingleGasNasaData(;name=nothing,
                       MM=nothing,
                       Hf=nothing,
                       H0=nothing,
                       Tlimit=nothing,
                       alow=nothing,
                       blow=nothing,
                       ahigh=nothing,
                       bhigh=nothing,
                       R=nothing) = 
                       new(name, MM,Hf,H0,Tlimit,alow,blow,ahigh,bhigh,R)
end


"""
    medium = SingleGasNasa(; mediumName     = Missing,
                             reference_p    = 101325,
                             reference_T    = 298.15,
                             p_default      = 101325,
                             T_default      = 293.15,
                             fluidConstants = nothing, 
                             fluidLimits    = FluidLimits(TMIN=200.0, TMAX=6000.0), 
                             data           = nothing)

Generate a `SingleGasNasa <: PureSubstance` medium object.
"""
struct SingleGasNasa <: PureSubstance
    infos::FluidInfos
    fluidConstants::SVector{1,IdealGasFluidConstants}
    fluidLimits::FluidLimits
    data::SingleGasNasaData

    function SingleGasNasa(; mediumName=nothing,
                             reference_p=101325,
                             reference_T=298.15,
                             p_default=101325,
                             T_default=293.15,
                             fluidConstants=nothing, 
                             fluidLimits=FluidLimits(TMIN=200.0, TMAX=6000.0), 
                             data=nothing)

        infos = FluidInfos(mediumName           = mediumName,
                           substanceNames       = [mediumName],
                           extraPropertiesNames = fill("",0),
                           ThermoStates         = IndependentVariables_pT,
                           baseProperties       = :BaseProperties_SingleGasNasa,
                           singleState          = false,
                           reducedX             = true,
                           fixedX               = false,
                           reference_p          = reference_p,
                           reference_T          = reference_T,
                           reference_X          = fill(1.0,1),
                           p_default            = p_default,
                           T_default            = T_default,
                           h_default            = h_T(data, T_default))

        fluidConstants.molarMass = data.MM
        fluidLimits.HMIN         = h_T(data, fluidLimits.TMIN)
        fluidLimits.HMAX         = h_T(data, fluidLimits.TMAX)

        new(infos, SVector{1,IdealGasFluidConstants}(fluidConstants), fluidLimits, data)
    end
end




### Functions specific for SingleGasNasa

"cp = cp_T(data::SingleGasNasaData, T) - Compute specific heat capacity at constant pressure from temperature and gas data"
cp_T(data::SingleGasNasaData, T::Float64)::Float64 = if T<data_Tlimit; data_R*
                                                (1/(T*T)*
                                                (data.alow[1]+T*
                                                (data.alow[2]+T*
                                                (1.0*data.alow[3]+T*
                                                (data.alow[4]+T*
                                                (data.alow[5]+T*(data.alow[6]+data.alow[7]*T))))))) else data.R*
                                                (1/(T*T)*
                                                (data.ahigh[1]+T*
                                                (data.ahigh[2]+T*
                                                (1.0*data.ahigh[3]+T*
                                                (data.ahigh[4]+T*
                                                (data.ahigh[5]+T*(data.ahigh[6]+data.ahigh[7]*T))))))) end

"""
    h = h_T(data, T, exclEnthForm=true, refChoice=ReferenceEnthalpy_ZeroAt0K, h_off=0.0)

Return specific enthalpy from temperature and gas data.

# Arguments
 
- `data::SingleGasNasaData`: Data of the SingleGasNasa medium.
- `T::Float64`: Temperature in [K].
- `exclEnthForm::Bool`: If true, enthalpy of formation Hf is not included in specific enthalpy h.
- `refChoice::ReferenceEnthalpy`: Choice of reference enthalpy.
- `h_off::Float64`: User defined offset for reference enthalpy, if referenceEnthalpy = ReferenceEnthalpy_UserDefined
"""
h_T(data::SingleGasNasaData, T::Float64, exclEnthForm::Bool=excludeEnthalpyOfFormation,
    refChoice::ReferenceEnthalpy=referenceEnthalpy, h_off::Float64=h_offset)::Float64 = 
                 (if T<data.Tlimit; data.R*
                 (
                 (-data.alow[1]+T*
                 (data.blow[1]+data.alow[2]*log(T)+T*
                 (1.0*data.alow[3]+T*
                 (0.5*data.alow[4]+T*
                 (1/3*data.alow[5]+T*(0.25*data.alow[6]+0.2*data.alow[7]*T))))))/T) else data.R*
                 (
                 (-data.ahigh[1]+T*
                 (data.bhigh[1]+data.ahigh[2]*log(T)+T*
                 (1.0*data.ahigh[3]+T*
                 (0.5*data.ahigh[4]+T*
                 (1/3*data.ahigh[5]+T*(0.25*data.ahigh[6]+0.2*data.ahigh[7]*T))))))/T) end)+(if exclEnthForm; -data.Hf else 0.0 end)+
                 (if refChoice==ReferenceEnthalpy_ZeroAt0K   ; data.H0 else 0.0 end)+
                 (if refChoice==ReferenceEnthalpy_UserDefined; h_off   else 0.0 end)


"Compute specific entropy from temperature and gas data"
s0_T(data::SingleGasNasaData,   # Ideal gas data
     T::Float64  # Temperature
    )::Float64 = (if T<data.Tlimit; data.R*
                 (data.blow[2]-0.5*data.alow[1]/(T*T)-data.alow[2]/T+data.alow[3]*log(T)+T*
                 (data.alow[4]+T*
                 (0.5*data.alow[5]+T*(1/3*data.alow[6]+0.25*data.alow[7]*T)))) else data.R*
                 (data.bhigh[2]-0.5*data.ahigh[1]/(T*T)-data.ahigh[2]/T+data.ahigh[3]*log(T)+T*
                 (data.ahigh[4]+T*
                 (0.5*data.ahigh[5]+T*(1/3*data.ahigh[6]+0.25*data.ahigh[7]*T)))) end)




### Set states
setState_pTX(m::SingleGasNasa,p,T,X) = ThermodynamicState_pT(p,T)
#setState_phX(m::SingleGasNasa,p,h,X) = ???     solve nonlinear equation with Brent algorithm (needs to be transferred from Modelica to Modia and stord in ModiaMath)
#setState_psX(m::SingleGasNasa,p,s,X) = ???     solve nonlinear equation with Brent algorithm
setState_dTX(m::SingleGasNasa,d,T,X) = ThermodynamicState_pT(d*m.data.R*T, T)

pressure(              m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = state.p
temperature(           m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = state.T
density(               m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = state.p/(m.data.R*state.T)
specificEnthalpy(      m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = h_T(m.data,state.T)
specificInternalEnergy(m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = h_T(m.data,state.T) - m.data.R*state.T
specificHeatCapacityCp(m::SingleGasNasa, state::ThermodynamicState_pT)::Float64 = s0_T(m.data,state.T) - m.data.R*log(state.p/m.infos.reference_p)


function standardCharacteristics(m::SingleGasNasa)::Dict{AbstractString,Any}
    p_ref = m.infos.reference_p
    T     = collect( range(m.fluidLimits.TMIN, stop=min(1600.0, m.fluidLimits.TMAX), length=501) )
    p     = [0.5e5, 1.0e5, 2.0e5]
    nT    = length(T)
    np    = length(p)
    h     = zeros(nT)
    u     = zeros(nT)
    cp    = zeros(nT)
    d     = zeros(nT,np)

    for i in 1:nT
        state = setState_pT(m,p_ref,T[i])
        h[i]  = specificEnthalpy(      m, state)
        u[i]  = specificInternalEnergy(m, state)
        cp[i] = specificHeatCapacityCp(m, state)
    end

    for j in 1:np
        for i in 1:nT
            d[i,j] = to_DensityDisplayUnit( density(m, ThermodynamicState_pT(p[j],T[i]) ) )     
        end
    end       

    mediumDict = Dict{AbstractString,Any}()
    mediumDict["T"]  = uconvert.(u"°C", T*1U"K")
    mediumDict["h"]  = h*1U"J/kg"
    mediumDict["u"]  = u*1U"J/kg"
    mediumDict["cp"] = cp*1U"J/(kg*K)"
    mediumDict["d(p=0.5 bar)"] = d[:,1]*1U"g/cm^3"
    mediumDict["d(p=1.0 bar)"] = d[:,2]*1U"g/cm^3"
    mediumDict["d(p=2.0 bar)"] = d[:,3]*1U"g/cm^3"
    return mediumDict
end


function standardPlot(m::SingleGasNasa; figure=1) 
    mediumDict = standardCharacteristics(m)
    ModiaMath.plot(mediumDict, [("h", "u"), "cp", ("d(p=0.5 bar)" ,
                                                   "d(p=1.0 bar)" ,
                                                   "d(p=2.0 bar)")], xAxis="T", heading=m.infos.mediumName, figure=figure)
end