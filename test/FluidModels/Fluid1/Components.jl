# License for this file: MIT
# Copyright 2018-2019, DLR Institute of System Dynamics and Control
#
"""
   module Fluid1_Components

Date: Jan. 22, 2019

Fluid component models with ModiaMedia based on the new approach described in 

- Zimmer, Bender, Pollok (2018): Robust Modeling of Directed Thermofluid Flows in Complex Networks.
  2nd Japanese Modelica Conference, Tokyo, Proceedings, pp 39-48. Download:
  https://www.modelica.org/events/modelica2018japan/conference-proceedings/modelica-final-proceedings-2018-Japan.pdf


# Main Authors of this module
```

      /|      Dirk Zimmer, Martin Otter
  ___/_|___   German Aerospace Center (DLR)
 /  /  /  /   Institute of System Dynamics and Control
/__/__/__/    https://www.dlr.de/sr/en/
   | /
   |/   
        
```
"""
module Fluid1_Components

using  Modia
using  ModiaMedia

using  ModiaMedia.StaticArrays # included via ModiaMedia, to avoid requirement to add it in the standard environment
using  ModiaMedia.Unitful      # included via ModiaMedia, to avoid requirement to add it in the standard environment


export add!
export BaseProperties
export StaticPressure, InertialPressure, SpecificEnthalpy, Temperature, MassFlowRate
export InPlug, OutPlug, FixedSource_pT, FixedSink, ShortPipe, Splitter, Junction, Compressor, ClosedVolume

# Variable types
StaticPressure(  ; args...) = Variable(; start=1.0e5, nominal=1.0e5, size = (), T = u"Pa"  , info = "Pressure for static mass-flow", args...)
InertialPressure(; args...) = Variable(; start=0.0  , nominal=1.0e5, size = (), T = u"Pa"  , info = "Pressure for accelerating the fluid", args...)
SpecificEnthalpy(; args...) = Variable(; start=1.0e4, nominal=1.0e4, size = (), T = u"J/kg", info = "Specific enthalpy", args...)
Temperature(     ; args...) = Variable(; start=300.0, nominal=300.0, size = (), T = u"K"   , info = "Temperature", args...)
MassFlowRate(    ; args...) = Variable(; start=0.0  ,                size = (), T = u"kg/s", info = "Mass flow rate", args...)


const MediumVariable() = Var(size=())
const MediumState()    = Var(size=())


const mediumFunctionMap = Dict{String, Symbol}("T" => :temperature,
                                               "d" => :density)

"""
    add!(result, portName, var::Vector{String})

Add the variables of `var` to port `portName` and store the vectors in `result`.
"""
function add!(result, portName::AbstractString, var::Vector{String})
    Medium = result[portName * ".Medium"][1]
    p = result[portName * ".p"]
    h = result[portName * ".h"]

    for j in eachindex(var)
        fc = getfield(ModiaMedia, mediumFunctionMap[var[j]])
        v  = fill(0.0, length(p))
        for i in eachindex(p)
            v[i] = fc(Medium, setState_ph(Medium,p[i],h[i]))
        end
        result[portName * "." * var[j]] = v
    end
 
    return nothing
end



###### Define BaseProperties models supported here (user can provide further definitions)

@model SimpleMedium_BaseProperties begin
    p_start = 1.0e5                # parameter
    T_start = 300.0                # parameter

    Medium = MediumVariable()      # input
    p      = StaticPressure(  start=p_start)  # input
    h      = SpecificEnthalpy(start=1.0)      # input
    T      = Temperature(     start=T_start)  # output
    d      = Float(start=1.0)      # output
    u      = Float(start=1.0)      # output
    der_d  = Float(start=1.0)      # output
    der_u  = Float(start=1.0)      # output
    state  = MediumState()         # output

@equations begin
    d     = density(Medium)
    u     = specificInternalEnergy_T(Medium,T)
    state = setState_pT(Medium,p,T)
    h     = specificEnthalpy(Medium,state)
    der_d = 0.0
    der_u = specificInternalEnergy_T_der_2(Medium,T)*der(T)
    end   
end


@model SimpleIdealGasMedium_BaseProperties begin
    p_start = 1.0e5                # parameter
    T_start = 300.0                # parameter

    Medium = MediumVariable()      # input
    p      = StaticPressure(  start=p_start)  # input
    h      = SpecificEnthalpy(start=1.0)      # input
    T      = Temperature(     start=T_start)  # output
    d      = Float(start=1.0)      # output
    u      = Float(start=1.0)      # output
    der_d  = Float(start=1.0)      # output
    der_u  = Float(start=1.0)      # output
    state  = MediumState()         # output

@equations begin
    d     = density_pT(Medium,p,T)
    u     = specificInternalEnergy_T(Medium,T)
    h     = specificEnthalpy_T(Medium,T)
    der_d = density_pT_der_2(Medium,p,T)*der(p) + density_pT_der_3(Medium,p,T)*der(T)
    der_u = specificInternalEnergy_T_der_2(Medium,T)*der(T)
    state = setState_pT(Medium,p,T)
    end   
end

BaseProperties(Medium::ModiaMedia.SimpleMedium        ; p_start=1e5, T_start=300.0, args...) = SimpleMedium_BaseProperties(        Medium=Medium, p_start=p_start, T_start=T_start, args...)
BaseProperties(Medium::ModiaMedia.SimpleIdealGasMedium; p_start=1e5, T_start=300.0, args...) = SimpleIdealGasMedium_BaseProperties(Medium=Medium, p_start=p_start, T_start=T_start, args...)


####### Define Fluid components

"""
    connector InPlug - Uni-directional flow of fluid (mass flow is entering the component)
"""
@model InPlug(:connector) begin
    Medium = MediumVariable()
    m_flow = MassFlowRate(flow = true)
    r      = InertialPressure()
    p      = StaticPressure() 
    h      = SpecificEnthalpy()
end




"""
    connector OutPlug - Uni-directional flow of fluid (mass flow is leaving the component)
"""
@model OutPlug(:connector) begin
    Medium = MediumVariable()
    m_flow = MassFlowRate(flow = true)
    r      = InertialPressure()
    p      = StaticPressure()
    h      = SpecificEnthalpy()
end



"""
    model FixedSource_pT - Fixed pressure and fixed specific enthalpy at source
"""
@model FixedSource_pT begin
    Medium  = MediumVariable()    # Medium MUST be redefined when instanting Source, otherwise error
    outPlug = OutPlug()
    p0 = StaticPressure(variability=parameter, info = "Fixed static pressure at source")
    T0 = Temperature(   variability=parameter, info = "Fixed temperature at source")
    @equations begin
        outPlug.Medium = Medium
        outPlug.r = 0.0
        outPlug.p = p0
        outPlug.h = specificEnthalpy(Medium, setState_pT(Medium,p0,T0))
    end
end



"""
    model FixedSink - Fixed pressure at sink
"""
@model FixedSink begin
    inPlug = InPlug()
    p0 = StaticPressure(variability=parameter, info = "Fixed static pressure at sink")
    @equations begin
        inPlug.p + inPlug.r = p0
    end
end



"""
    model ShortPipe - Simple pipe (without volume)
"""
@model ShortPipe begin
    inPlug  = InPlug()
    outPlug = OutPlug()

    k = Parameter(1.0e5, info="Quadratic pressure drop coefficient")
    l = Parameter(1.0  , info="Length of pipe")
    A = Parameter(1.0  , info="Area of pipe")    # should be taken from a global world model or propagated through the connection structure

    dp     = StaticPressure(  info = "delta pressure in flow direction")
    dr     = InertialPressure(info = "delta inertial pressure in flow direction")
    dh     = SpecificEnthalpy(info = "delta in spec. enthalpy in flow direction")
    m_flow = MassFlowRate()

    state = Var(size=())
    eta   = Float(size=())
    @equations begin
        # Medium propagation
        outPlug.Medium = inPlug.Medium

        # mass flow balance
        m_flow = inPlug.m_flow
        inPlug.m_flow + outPlug.m_flow = 0

        # Propagation of specific quantities
        outPlug.p  = inPlug.p + dp;
        outPlug.r  = inPlug.r + dr;
        outPlug.h  = inPlug.h + dh;

        dp = -m_flow*abs(m_flow) * k    # pressure drop
        dh = 0.0                        # isenthalpic process

        # inertial component of pressure
        der(m_flow)/A*l = -dr

        # compute dynamic viscosity to test dispatch of function
        state = setState_ph(inPlug.Medium, inPlug.p, inPlug.h)
        eta   = dynamicViscosity(inPlug.Medium, state)
    end
end



"""
    model Splitter - one inflow, two outflows
"""
@model Splitter begin
    inPlugA  = InPlug()
    outPlugB = OutPlug()
    outPlugC = OutPlug()

    eps = 0.000001

    @equations begin
        # Medium propagation
        outPlugB.Medium = inPlugA.Medium
        outPlugC.Medium = inPlugA.Medium

        # mass flow balance
        inPlugA.m_flow + outPlugB.m_flow + outPlugC.m_flow = 0.0

        # propagation of flow and pressure balance
        outPlugB.r = inPlugA.r
        outPlugB.p = inPlugA.p
        outPlugB.h = inPlugA.h

        outPlugC.r = inPlugA.r
        outPlugC.p = inPlugA.p
        outPlugC.h = inPlugA.h
    end
end


"""
    model Junction - two inflows, one outflow
"""
@model Junction begin
    inPlugA  = InPlug()
    inPlugB  = InPlug()
    outPlugC = OutPlug()

    eps = 0.000001

    @equations begin
        # Medium propagation
        inPlugA.Medium = outPlugC.Medium
        inPlugB.Medium = outPlugC.Medium

        # mass flow balance
        inPlugA.m_flow + inPlugB.m_flow + outPlugC.m_flow = 0.0

        # total pressure balance
        inPlugA.p + inPlugA.r = outPlugC.p + outPlugC.r
        inPlugB.p + inPlugB.r = outPlugC.p + outPlugC.r


        # mixing of fluid
            # formulation safeguarded against zero-massflow and backward flow
            (abs(outPlugC.m_flow) + 2*eps)*outPlugC.h = (abs(inPlugA.m_flow) + eps)*inPlugA.h + (abs(inPlugB.m_flow) + eps)*inPlugB.h

            # for simplicity m_flow is taken for the weighted mean of p but it should be V_flow
            (abs(outPlugC.m_flow) + 2*eps)*outPlugC.p = (abs(inPlugA.m_flow) + eps)*inPlugA.p + (abs(inPlugB.m_flow) + eps)*inPlugB.p
    end
end



"""
    model Compressor - Simple compressor
"""
@model Compressor begin
    inPlug  = InPlug()
    outPlug = OutPlug()

    pressureRatio = Parameter(2.0, info="Fixed pressure ratio")
    efficiency    = Parameter(0.9)
    kappa         = Parameter(1.4)
    R             = Parameter(8.314, info="Gas constant")
    
    k = Parameter(1.0e5, info="Quadratic pressure drop coefficient")
    l = Parameter(1.0  , info="Length of compressor")
    A = Parameter(1.0  , info="Area of compressor")    # should be taken from a global world model or propagated through the connection structure

    dr     = InertialPressure(info = "delta inertial pressure in flow direction")
    dh     = SpecificEnthalpy(info = "delta in spec. enthalpy in flow direction")
    m_flow = MassFlowRate()
    y      = Float(info="Specific energy inserted into the mass-flow")
    k_frac = Float()

    @equations begin
        # Medium propagation
        outPlug.Medium = inPlug.Medium
        
        k_frac = kappa/(kappa-1.0)

        # mass flow balance
        m_flow = inPlug.m_flow
        inPlug.m_flow + outPlug.m_flow = 0

        # Propagation of specific quantities
        outPlug.p = pressureRatio*inPlug.p
        outPlug.r = inPlug.r + dr
        outPlug.h = inPlug.h + dh
        y  = k_frac*R*(pressureRatio^(1/k_frac)-1)
        dh = y/efficiency;

        # inertial component of pressure
        der(m_flow)/A*l = -dr
    end
end



"""
    model ClosedVolume - Simple volume
"""
@model ClosedVolume begin
    Medium  = AbstractMedium    # Medium MUST be redefined when instanting ClosedVolume, otherwise error

    inPlug  = InPlug()
    outPlug = OutPlug()

    V  = 1.0       # Fixed volume
    p0 = Medium.infos.p_default     # initial pressure
    T0 = Medium.infos.T_default     # initial temperature

    M     = Float(info="Mass of volume"      , start=1.0)
    der_M = Float(info="Total mass flow rate", start=0.0)

    #medium = getBaseProperties(Medium)(p_start=p0, T_start=T0)
    medium = BaseProperties(Medium; p_start=p0, T_start=T0)

    h = Float(size=())
@equations begin
    # Medium propagation
    # outPlug.Medium = inPlug.Medium
    outPlug.Medium = Medium


    # Mass balance
    der_M          = inPlug.m_flow + outPlug.m_flow
    medium.der_d*V = der_M

    # Energy balance
    M = medium.d*V
    der_M*medium.u + M*medium.der_u = inPlug.m_flow*(inPlug.h + inPlug.r/medium.d) + outPlug.m_flow*outPlug.h

    # Properties at outPlug
    outPlug.p = medium.p
    outPlug.h = medium.h
    outPlug.r = 0.0

    # Properties of inPlug
    inPlug.p + inPlug.r = medium.p

    # For testing 
    h = specificEnthalpy(Medium, medium.state)
    end
end


end