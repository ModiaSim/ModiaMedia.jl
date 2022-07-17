#
# This file is part of module 
#   ModiaMedia (ModiaMedia/src/ModiaMedia.jl)
#

# Dictionary of medium data


const SingleGasNasa_excludeEnthalpyOfFormation = true                         # If true, enthalpy of formation Hf is not included in specific enthalpy hc
const SingleGasNasa_referenceEnthalpy          = ReferenceEnthalpy_ZeroAt0K   # Choice of reference enthalpy
const SingleGasNasa_h_offset                   = 0.0                          # User defined offset for reference enthalpy, if ReferenceEnthalpy = ReferenceEnthalpy_UserDefined




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
end


"""
    state = SingleGasNasaState(Medium, p, T)

Generate a `SingleGasNasaState <: ThermodynamicState` object containing
pressure `p` [Pa] and temperature `T` [K] as thermodynamic states.
"""
mutable struct SingleGasNasaState <: ThermodynamicState
    Medium::SingleGasNasa
    p::Float64
    T::Float64
end


function SingleGasNasa(; mediumName=nothing,
                         reference_p=101325.0,
                         reference_T=298.15,
                         p_default=101325.0,
                         T_default=293.15,
                         fluidConstants=nothing, 
                         fluidLimits=FluidLimits(TMIN=200.0, TMAX=6000.0), 
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

    fluidConstants.molarMass = data.MM
    infos.h_default          = h_T(data, T_default)
    fluidLimits.HMIN         = h_T(data, fluidLimits.TMIN)
    fluidLimits.HMAX         = h_T(data, fluidLimits.TMAX)

    Medium = SingleGasNasa(infos, SVector{1,IdealGasFluidConstants}(fluidConstants), fluidLimits, data)

    return Medium
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
- `h_off::Float64`: User defined offset for reference enthalpy, if SingleGasNasa_referenceEnthalpy = ReferenceEnthalpy_UserDefined
"""
h_T(data::SingleGasNasaData, T::Float64; 
    exclEnthForm::Bool=SingleGasNasa_excludeEnthalpyOfFormation,
    refChoice::ReferenceEnthalpy=SingleGasNasa_referenceEnthalpy, 
    h_off::Float64=SingleGasNasa_h_offset)::Float64 = 
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
                 (if refChoice==ReferenceEnthalpy_ZeroAt0K   ; data.H0 else 0.0 end) +
                 (if refChoice==ReferenceEnthalpy_UserDefined; h_off   else 0.0 end)


"Compute specific enthalpy, low T region; reference is decided by the refChoice input, or by the referenceChoice package constant by default"
h_Tlow(data::SingleGasNasaData,      # Ideal gas data
    T::Float64;                      # Temperature
    exclEnthForm::Bool=SingleGasNasa_excludeEnthalpyOfFormation,
    refChoice::ReferenceEnthalpy=SingleGasNasa_referenceEnthalpy, 
    h_off::Float64=SingleGasNasa_h_offset) = 
        data.R*((-data.alow[1]+T*(data.blow[1]+data.alow[2]*log(T)+T*(1.0*data.alow[3]+T*(0.5*data.alow[4]+T*(1/3*data.alow[5]+T*(0.25*data.alow[6]+0.2*data.alow[7]*T))))))/T)+
             (if exclEnthForm; -data.Hf else 0.0 end)+
             (if refChoice==ReferenceEnthalpy_ZeroAt0K   ; data.H0 else 0.0 end) +
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



"Dynamic viscosity of low pressure gases"
function dynamicViscosityLowPressure(
    T::Float64,                      # Gas temperature
    Tc::Float64,                     # Critical temperature of gas
    M::Float64,                      # Molar mass of gas
    Vc::Float64,                     # Critical molar volume of gas
    w::Float64,                      # Acentric factor of gas
    mu::Float64;                     # Modelica.Media.Interfaces.Types.DipoleMoment; Dipole moment of gas molecule
    k::Float64=0.0,                  # Special correction for highly polar substances
    # returns eta::Float64,          # Dynamic viscosity of gas
    )
    Const1_SI::Float64 = 40.785*10.0^(-9.5)  # Constant in formula for eta converted to SI units
    Const2_SI::Float64 = 131.3/1000.0        # Constant in formula for mur converted to SI units
    mur::Float64       = Const2_SI*mu/sqrt(Vc*Tc) # Dimensionless dipole moment of gas molecule
    Fc::Float64        = 1 - 0.2756*w + 0.059035*mur^4 + k                 # Factor to account for molecular shape and polarities of gas
    Tstar = 1.2593*T/Tc
    Ov = 1.16145*Tstar^(-0.14874)+0.52487*exp(-0.7732*Tstar)+2.16178*exp(-2.43787*Tstar)
    eta = Const1_SI*Fc*sqrt(M*T)/(Vc^(2/3)*Ov)
    return eta
end

T_h( data::SingleGasNasaData, h::Float64) = ModiaMedia.solveOneNonlinearEquation(T->h-h_T(data,T), 200.0, 6000.0; u_nominal=300.0)
T_ps(m::SingleGasNasa, p::Float64, s::Float64) = ModiaMedia.solveOneNonlinearEquation(T->s0_T(m.data,T)-m.data.R*log(p/m.infos.reference_p), 200.0, 6000.0; u_nominal=300.0)


### Set states
setState_pTX(m::SingleGasNasa,p,T,X) = SingleGasNasaState(m,p,T)
setState_phX(m::SingleGasNasa,p,h,X) = SingleGasNasaState(m,p,T_h(m.data,h))
setState_psX(m::SingleGasNasa,p,s,X) = SingleGasNasaState(m,p,T_ps(m,p,s))
setState_dTX(m::SingleGasNasa,d,T,X) = SingleGasNasaState(m,d*m.data.R*T, T)
isenthalpicState(m::SingleGasNasa, state::SingleGasNasaState, dp::Float64) = SingleGasNasaState(m, state.p+dp, state.T)


setState_pTX!(state::SingleGasNasaState,p,T,X) = begin state.p=p; state.T=T; nothing end
setState_phX!(state::SingleGasNasaState,p,h,X) = begin state.p=p; state.T=T_h(state.Medium.data,h); nothing end
setState_psX!(state::SingleGasNasaState,p,s,X) = begin state.p=p; state.T=T_ps(state.Medium,p,s); nothing end
setState_dTX!(state::SingleGasNasaState,d,T,X) = begin state.p=d*state.Medium.data.R*T; state.T=T; nothing end
isenthalpicState!(state_b::SingleGasNasaState, state_a::SingleGasNasaState, dp::Float64) = begin state_b.p = state_a.p+dp; state_b.T = state_a.T; nothing end


pressure(              m::SingleGasNasa, state::SingleGasNasaState)::Float64 = state.p
temperature(           m::SingleGasNasa, state::SingleGasNasaState)::Float64 = state.T
density(               m::SingleGasNasa, state::SingleGasNasaState)::Float64 = state.p/(m.data.R*state.T)
specificEnthalpy(      m::SingleGasNasa, state::SingleGasNasaState)::Float64 = h_T(m.data,state.T)
specificInternalEnergy(m::SingleGasNasa, state::SingleGasNasaState)::Float64 = h_T(m.data,state.T) - m.data.R*state.T
specificHeatCapacityCp(m::SingleGasNasa, state::SingleGasNasaState)::Float64 = s0_T(m.data,state.T) - m.data.R*log(state.p/m.infos.reference_p)

function dynamicViscosity(m::SingleGasNasa, state::SingleGasNasaState)::Float64
    @assert(m.fluidConstants[1].hasCriticalData, "Failed to compute dynamicViscosity: For the species \""+m.mediumName+"\" no critical data is available.")
    @assert(m.fluidConstants[1].hasDipoleMoment, "Failed to compute dynamicViscosity: For the species \""+m.mediumName+"\" no critical data is available.")
    c   = m.fluidConstants[1]
    eta = dynamicViscosityLowPressure(state.T, c.criticalTemperature, 
                                               c.molarMass, 
                                               c.criticalMolarVolume, 
                                               c.acentricFactor, 
                                               c.dipoleMoment)
    return eta
end


function standardCharacteristics(m::SingleGasNasa)
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
            d[i,j] = to_DensityDisplayUnit( density(SingleGasNasaState(m,p[j],T[i]) ) )     
        end
    end       

#=
    mediumDict = Dict{AbstractString,Any}()
    mediumDict["T"]  = uconvert.(u"°C", T*1u"K")
    mediumDict["h"]  = h*1u"J/kg"
    mediumDict["u"]  = u*1u"J/kg"
    mediumDict["cp"] = cp*1u"J/(kg*K)"
    mediumDict["d(p=0.5 bar)"] = d[:,1]*1u"g/cm^3"
    mediumDict["d(p=1.0 bar)"] = d[:,2]*1u"g/cm^3"
    mediumDict["d(p=2.0 bar)"] = d[:,3]*1u"g/cm^3"
=#
    
    mediumSignalTable = SignalTable(
        "T"  => Var(values = ustrip.(uconvert.(u"°C", T*1u"K")), unit="°C", independent=true),
        "h"  => Var(values =  h, unit ="J/kg"),
        "u"  => Var(values =  u, unit ="J/kg"),
        "cp" => Var(values = cp, unit = "J/(kg*K)")        
    )    
    return mediumSignalTable
end


function standardPlot(m::SingleGasNasa, plotFunction::Function; figure=1) 
    mediumSignalTable = standardCharacteristics(m)
    #plot(mediumDict, [("h", "u"), "cp", ("d(p=0.5 bar)" ,
    #                  "d(p=1.0 bar)" ,
    #                  "d(p=2.0 bar)")], xAxis="T", heading=m.infos.mediumName, figure=figure)
    plotFunction(mediumSignalTable, [("h", "u"), "cp"], xAxis="T", heading=m.infos.mediumName, figure=figure)    
end