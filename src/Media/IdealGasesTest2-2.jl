
"Data and models of ideal gases (single, fixed and dynamic mixtures) from NASA source"
module IdealGases

    "Common packages and data for the ideal gas models"
    module Common

        "Coefficient data record for properties of ideal gases based on NASA source"
        mutable struct DataRecord
            name::String                     # Name of ideal gas
            MM::Float64                      # Molar mass
            Hf::Float64                      # Enthalpy of formation at 298.15K
            H0::Float64                      # H0(298.15K) - H0(0K)
            Tlimit::Float64                  # Temperature limit between low and high data sets
            alow::Array{Float64, 1}          # Low temperature coefficients a
            blow::Array{Float64, 1}          # Low temperature constants b
            ahigh::Array{Float64, 1}         # High temperature coefficients a
            bhigh::Array{Float64, 1}         # High temperature constants b
            R::Float64                       # Gas constant
        end

        "Basic Functions for ideal gases: cp, h, s, thermal conductivity, viscosity"
        module Functions
            const excludeEnthalpyOfFormation=true   # If true, enthalpy of formation Hf is not included in specific enthalpy h
            const referenceChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K   # Choice of reference enthalpy
            const h_offset=0.0   # User defined offset for reference enthalpy, if referenceChoice = UserDefined
            const methodForThermalConductivity = 1

            "Compute specific heat capacity at constant pressure from temperature and gas data"
            function cp_T(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                # returns cp::Float64,           # Specific heat capacity at temperature T
                )
                cp = smooth(0, if T<data.Tlimit; data.R*
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
                    (data.ahigh[5]+T*(data.ahigh[6]+data.ahigh[7]*T))))))) end)
                
                return cp
            end

            "Compute specific heat capacity at constant pressure, low T region"
            function cp_Tlow(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                # returns cp::Float64,           # Specific heat capacity at temperature T
                )
                cp = data.R*
                    (1/(T*T)*
                    (data.alow[1]+T*
                    (data.alow[2]+T*
                    (1.0*data.alow[3]+T*
                    (data.alow[4]+T*
                    (data.alow[5]+T*(data.alow[6]+data.alow[7]*T)))))))
                
                return cp
            end

            "Compute specific heat capacity at constant pressure, low T region"
            function cp_Tlow_der(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                dT::Float64,                     # Temperature derivative
                # returns cp_der::Float64,       # Derivative of specific heat capacity
                )
                cp_der = dT*data.R/(T*T*T)*
                    (-2*data.alow[1]+T*
                    (-data.alow[2]+T*T*
                    (data.alow[4]+T*
                    (2.0*data.alow[5]+T*(3.0*data.alow[6]+4.0*data.alow[7]*T)))))
                
                return cp_der
            end

            "Compute specific enthalpy from temperature and gas data; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            function h_T(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                exclEnthForm::Bool=excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::Modelica.Media.Interfaces.Choices.ReferenceEnthalpy=referenceChoice,   # Choice of reference enthalpy
                h_off::Float64=h_offset,         # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                # returns h::Float64,            # Specific enthalpy at temperature T
                )
                h = smooth(0, 
                    (if T<data.Tlimit; data.R*
                    (
                    (-data.alow[1]+T*
                    (data.blow[1]+data.alow[2]*Math.log(T)+T*
                    (1.0*data.alow[3]+T*
                    (0.5*data.alow[4]+T*
                    (1/3*data.alow[5]+T*(0.25*data.alow[6]+0.2*data.alow[7]*T))))))/T) else data.R*
                    (
                    (-data.ahigh[1]+T*
                    (data.bhigh[1]+data.ahigh[2]*Math.log(T)+T*
                    (1.0*data.ahigh[3]+T*
                    (0.5*data.ahigh[4]+T*
                    (1/3*data.ahigh[5]+T*(0.25*data.ahigh[6]+0.2*data.ahigh[7]*T))))))/T) end)+(if exclEnthForm; -data.Hf else 0.0 end)+
                    (if     (refChoice==Choices.ReferenceEnthalpy.ZeroAt0K); data.H0 else 0.0 end)+
                    (if refChoice==Choices.ReferenceEnthalpy.UserDefined; h_off else 0.0 end))
                
                return h
            end

            "Derivative function for h_T"
            function h_T_der(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                exclEnthForm::Bool=excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::Modelica.Media.Interfaces.Choices.ReferenceEnthalpy=referenceChoice,   # Choice of reference enthalpy
                h_off::Float64=h_offset,         # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                dT::Float64,                     # Temperature derivative
                # returns h_der::Float64,        # Specific enthalpy at temperature T
                )
                h_der = dT*Modelica.Media.IdealGases.Common.Functions.cp_T(data, T)
                
                return h_der
            end

            "Compute specific enthalpy, low T region; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            function h_Tlow(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                exclEnthForm::Bool=excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::Modelica.Media.Interfaces.Choices.ReferenceEnthalpy=referenceChoice,   # Choice of reference enthalpy
                h_off::Float64=h_offset,         # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                # returns h::Float64,            # Specific enthalpy at temperature T
                )
                h = data.R*
                    (
                    (-data.alow[1]+T*
                    (data.blow[1]+data.alow[2]*Math.log(T)+T*
                    (1.0*data.alow[3]+T*
                    (0.5*data.alow[4]+T*
                    (1/3*data.alow[5]+T*(0.25*data.alow[6]+0.2*data.alow[7]*T))))))/T)+(if exclEnthForm; -data.Hf else 0.0 end)+
                    (if     (refChoice==Choices.ReferenceEnthalpy.ZeroAt0K); data.H0 else 0.0 end)+
                    (if refChoice==Choices.ReferenceEnthalpy.UserDefined; h_off else 0.0 end)
                
                return h
            end

            "Compute specific enthalpy, low T region; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            function h_Tlow_der(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                exclEnthForm::Bool=excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::Modelica.Media.Interfaces.Choices.ReferenceEnthalpy=referenceChoice,   # Choice of reference enthalpy
                h_off::Float64=h_offset,         # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                dT::Float64,                     # Temperature derivative
                # returns h_der::Float64,        # Derivative of specific enthalpy at temperature T
                )
                h_der = dT*Modelica.Media.IdealGases.Common.Functions.cp_Tlow(data, T)
                
                return h_der
            end

            "Compute specific entropy from temperature and gas data"
            function s0_T(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                # returns s::Float64,            # Specific entropy at temperature T
                )
                s = if T<data.Tlimit; data.R*
                    (data.blow[2]-0.5*data.alow[1]/(T*T)-data.alow[2]/T+data.alow[3]*Math.log(T)+T*
                    (data.alow[4]+T*
                    (0.5*data.alow[5]+T*(1/3*data.alow[6]+0.25*data.alow[7]*T)))) else data.R*
                    (data.bhigh[2]-0.5*data.ahigh[1]/(T*T)-data.ahigh[2]/T+data.ahigh[3]*Math.log(T)+T*
                    (data.ahigh[4]+T*
                    (0.5*data.ahigh[5]+T*(1/3*data.ahigh[6]+0.25*data.ahigh[7]*T)))) end
                
                return s
            end

            "Compute specific entropy, low T region"
            function s0_Tlow(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                # returns s::Float64,            # Specific entropy at temperature T
                )
                s = data.R*
                    (data.blow[2]-0.5*data.alow[1]/(T*T)-data.alow[2]/T+data.alow[3]*Math.log(T)+T*
                    (data.alow[4]+T*
                    (0.5*data.alow[5]+T*(1/3*data.alow[6]+0.25*data.alow[7]*T))))
                
                return s
            end

            "Compute derivative of specific entropy, low T region"
            function s0_Tlow_der(
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                T::Float64,                      # Temperature
                T_der::Float64,                  # Temperature derivative
                # returns s::Float64,            # Specific entropy at temperature T
                )
                s = data.R*
                    (data.blow[2]-0.5*data.alow[1]/(T*T)-data.alow[2]/T+data.alow[3]*Math.log(T)+T*
                    (data.alow[4]+T*
                    (0.5*data.alow[5]+T*(1/3*data.alow[6]+0.25*data.alow[7]*T))))
                
                return s
            end
#=
            "Dynamic viscosity of low pressure gases"
            function dynamicViscosityLowPressure(
                T::Float64,                      # Gas temperature
                Tc::Float64,                     # Critical temperature of gas
                M::Float64,                      # Molar mass of gas
                Vc::Float64,                     # Critical molar volume of gas
                w::Float64,                      # Acentric factor of gas
                mu::Modelica.Media.Interfaces.Types.DipoleMoment,   # Dipole moment of gas molecule
                k::Float64=0.0,                  # Special correction for highly polar substances
                # returns eta::Float64,          # Dynamic viscosity of gas
                Const1_SI::Float64=40.785*10^(-9.5),   # Constant in formula for eta converted to SI units
                Const2_SI::Float64=131.3/1000.0,   # Constant in formula for mur converted to SI units
                mur::Float64=Const2_SI*mu/sqrt(Vc*Tc),   # Dimensionless dipole moment of gas molecule
                Fc::Float64=1-0.2756*w+0.059035*mur^4+k,   # Factor to account for molecular shape and polarities of gas
                Tstar::Float64,                  # Dimensionless temperature defined by equation below
                Ov::Float64,                     # Viscosity collision integral for the gas
                )
                Tstar = 1.2593*T/Tc
                                Ov = 1.16145*Tstar^(-0.14874)+0.52487*Modelica.Math.exp(-0.7732*Tstar)+2.16178*Modelica.Math.exp(-2.43787*Tstar)
                                eta = Const1_SI*Fc*sqrt(M*T)/(Vc^(2/3)*Ov)
                
                return 
            end

            "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
            function thermalConductivityEstimate(
                Cp::Modelica.Media.Interfaces.Types.SpecificHeatCapacity,   # Constant pressure heat capacity
                eta::Modelica.Media.Interfaces.Types.DynamicViscosity,   # Dynamic viscosity
                method::Int64 = 1,               # 1: Eucken Method, 2: Modified Eucken Method
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                # returns lambda::Modelica.Media.Interfaces.Types.ThermalConductivity,   # Thermal conductivity [W/(m.k)]
                )
                lambda = if method==1; eta*(Cp-data.R+(9/4)*data.R) else eta*(Cp-data.R)*
                    (1.32+1.77/((Cp/Modelica.Constants.R)-1.0)) end
                
                return lambda
            end
        end
=#
        "Medium model of an ideal gas based on NASA source"
        module SingleGasNasa
            mutable struct ThermodynamicState
                p::AbsolutePressure              # Absolute pressure of medium
                T::Temperature                   # Temperature of medium
            end
#=
            const data::IdealGases.Common.DataRecord   # Data record of ideal gas substance
            const fluidConstants::Array{FluidConstants, 1}   # Constant data for the fluid
=#
            @model BaseProperties begin
                @inherits assert, T, String, mediumName, MM, data, R, h, Modelica, u, d, p, state
                @equations begin
                    assert(T>=200 && T<=6000, "
                    Temperature T (= "+String(T)+" K) is not in the allowed range
                    200 K <= T <= 6000 K required from medium model \""+mediumName+"\".
                    ")
                    MM = data.MM
                    R = data.R
                    h = Modelica.Media.IdealGases.Common.Functions.h_T(data, T, Modelica.Media.IdealGases.Common.Functions.excludeEnthalpyOfFormation, Modelica.Media.IdealGases.Common.Functions.referenceChoice, Modelica.Media.IdealGases.Common.Functions.h_offset)
                    u = h-R*T
                    d = p/(R*T)
                    state.T = T
                    state.p = p
                    end
            end

            "Return thermodynamic state as function of p, T and composition X"
            function setState_pTX(
                p::AbsolutePressure,             # Pressure
                T::Temperature,                  # Temperature
                X::Array{MassFraction, 1}=reference_X,   # Mass fractions
                # returns state::ThermodynamicState,
                )
                state = ThermodynamicState(p=p, T=T)
                
                return state
            end

            "Return thermodynamic state as function of p, h and composition X"
            function setState_phX(
                p::AbsolutePressure,             # Pressure
                h::SpecificEnthalpy,             # Specific enthalpy
                X::Array{MassFraction, 1}=reference_X,   # Mass fractions
                # returns state::ThermodynamicState,
                )
                state = ThermodynamicState(p=p, T=T_h(h))
                
                return state
            end

            "Return thermodynamic state as function of p, s and composition X"
            function setState_psX(
                p::AbsolutePressure,             # Pressure
                s::SpecificEntropy,              # Specific entropy
                X::Array{MassFraction, 1}=reference_X,   # Mass fractions
                # returns state::ThermodynamicState,
                )
                state = ThermodynamicState(p=p, T=T_ps(p, s))
                
                return state
            end

            "Return thermodynamic state as function of d, T and composition X"
            function setState_dTX(
                d::Density,                      # Density
                T::Temperature,                  # Temperature
                X::Array{MassFraction, 1}=reference_X,   # Mass fractions
                # returns state::ThermodynamicState,
                )
                state = ThermodynamicState(p=d*data.R*T, T=T)
                
                return state
            end
            #=
            function setSmoothState(
                )
                state = ThermodynamicState(p=Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T=Media.Common.smoothStep(x, state_a.T, state_b.T, x_small))
                
                return 
            function pressure(
                )
                p = state.p
                
                return 
            function temperature(
                )
                T = state.T
                
                return 
            function density(
                )
                d = state.p/(data.R*state.T)
                
                return 
            function specificEnthalpy(
                )
                h = Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T)
                
                return 
            function specificInternalEnergy(
                )
                u = Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T)-data.R*state.T
                
                return 
            function specificEntropy(
                )
                s = Modelica.Media.IdealGases.Common.Functions.s0_T(data, state.T)-data.R*Modelica.Math.log(state.p/reference_p)
                
                return 
            function specificGibbsEnergy(
                )
                g = Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T)-state.T*specificEntropy(state)
                
                return 
            function specificHelmholtzEnergy(
                )
                f = Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T)-data.R*state.T-state.T*specificEntropy(state)
                
                return 
            function specificHeatCapacityCp(
                )
                cp = Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T)
                
                return 
            function specificHeatCapacityCv(
                )
                cv = Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T)-data.R
                
                return 
            function isentropicExponent(
                )
                gamma = specificHeatCapacityCp(state)/specificHeatCapacityCv(state)
                
                return 
            function velocityOfSound(
                )
                a = sqrt(max(0, data.R*state.T*Modelica.Media.IdealGases.Common.Functions.cp_T(data, state.T)/specificHeatCapacityCv(state)))
                
                return 
=#
            "Approximate method of calculating h_is from upstream properties and downstream pressure"
            function isentropicEnthalpyApproximation(
                p2::Float64,                     # Downstream pressure
                state::ThermodynamicState,       # Properties at upstream location
                exclEnthForm::Bool=Functions.excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::ReferenceEnthalpy=Functions.referenceChoice,   # Choice of reference enthalpy
                h_off::SpecificEnthalpy=Functions.h_offset,   # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                # returns h_is::Float64,         # Isentropic enthalpy
                gamma::IsentropicExponent=isentropicExponent(state),   # Isentropic exponent
                )
                h_is = Modelica.Media.IdealGases.Common.Functions.h_T(data, state.T, exclEnthForm, refChoice, h_off)+gamma/(gamma-1.0)*state.p/density(state)*((p2/state.p)^((gamma-1)/gamma)-1.0)
                
                return 
            end
#=
            function isentropicEnthalpy(
                exclEnthForm::Bool=Functions.excludeEnthalpyOfFormation,   # If true, enthalpy of formation Hf is not included in specific enthalpy h
                refChoice::ReferenceEnthalpy=Functions.referenceChoice,   # Choice of reference enthalpy
                h_off::SpecificEnthalpy=Functions.h_offset,   # User defined offset for reference enthalpy, if referenceChoice = UserDefined
                )
                h_is = isentropicEnthalpyApproximation(p_downstream, refState, exclEnthForm, refChoice, h_off)
                
                return 
            function isobaricExpansionCoefficient(
                )
                beta = 1/state.T
                
                return 
            function isothermalCompressibility(
                )
                kappa = 1.0/state.p
                
                return 
            function density_derp_T(
                )
                ddpT = 1/(state.T*data.R)
                
                return 
            function density_derT_p(
                )
                ddTp = -state.p/(state.T*state.T*data.R)
                
                return 
            function density_derX(
                )
                dddX = fill(0, nX)
                
                return 
            function dynamicViscosity(
                )
                assert(fluidConstants[1].hasCriticalData, "Failed to compute dynamicViscosity: For the species \""+mediumName+"\" no critical data is available.")
                                assert(fluidConstants[1].hasDipoleMoment, "Failed to compute dynamicViscosity: For the species \""+mediumName+"\" no critical data is available.")
                                eta = Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(state.T, fluidConstants[1].criticalTemperature, fluidConstants[1].molarMass, fluidConstants[1].criticalMolarVolume, fluidConstants[1].acentricFactor, fluidConstants[1].dipoleMoment)
                
                return 
            function thermalConductivity(
                method::Int64=Functions.methodForThermalConductivity,   # 1: Eucken Method, 2: Modified Eucken Method
                )
                assert(fluidConstants[1].hasCriticalData, "Failed to compute thermalConductivity: For the species \""+mediumName+"\" no critical data is available.")
                                lambda = Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(specificHeatCapacityCp(state), dynamicViscosity(state), method=method, data=data)
                
                return 
            function molarMass(
                )
                MM = data.MM
                
                return 
=#
            "Compute temperature from specific enthalpy"
            function T_h(
                h::SpecificEnthalpy,             # Specific enthalpy
                # returns T::Temperature,        # Temperature
#=
                "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
                module Internal
                    mutable struct f_nonlinear_Data
                    function f_nonlinear(
                        )
                        y = Modelica.Media.IdealGases.Common.Functions.h_T(f_nonlinear_data, x)
                        
                        return 
                    function solve(
                        return 
                end
=#                )
                T = Internal.solve(h, 200, 6000, 1.0e5, [1], data)
                
                return 
            end

            "Compute temperature from pressure and specific entropy"
            function T_ps(
                p::AbsolutePressure,             # Pressure
                s::SpecificEntropy,              # Specific entropy
                # returns T::Temperature,        # Temperature

                "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
#=
                module Internal
                    mutable struct f_nonlinear_Data
                    function f_nonlinear(
                        )
                        y = Modelica.Media.IdealGases.Common.Functions.s0_T(f_nonlinear_data, x)-data.R*Modelica.Math.log(p/reference_p)
                        
                        return 
                    function solve(
                        return 
                end
=#
                )
                T = Internal.solve(s, 200, 6000, p, [1], data)
                
                return 
            end

            "Dynamic viscosity of low pressure gases"
            function dynamicViscosityLowPressure(
                T::Float64,                      # Gas temperature
                Tc::Float64,                     # Critical temperature of gas
                M::Float64,                      # Molar mass of gas
                Vc::Float64,                     # Critical molar volume of gas
                w::Float64,                      # Acentric factor of gas
                mu::Modelica.Media.Interfaces.Types.DipoleMoment,   # Dipole moment of gas molecule
                k::Float64=0.0,                  # Special correction for highly polar substances
                # returns eta::Modelica.Media.Interfaces.Types.DynamicViscosity,   # Dynamic viscosity of gas
                Const1_SI::Float64=40.785*10^(-9.5),   # Constant in formula for eta converted to SI units
                Const2_SI::Float64=131.3/1000.0,   # Constant in formula for mur converted to SI units
                mur::Float64=Const2_SI*mu/sqrt(Vc*Tc),   # Dimensionless dipole moment of gas molecule
                Fc::Float64=1-0.2756*w+0.059035*mur^4+k,   # Factor to account for molecular shape and polarities of gas
                Tstar::Float64,                  # Dimensionless temperature defined by equation below
                Ov::Float64,                     # Viscosity collision integral for the gas
                )
                eta = Functions.dynamicViscosityLowPressure(T, Tc, M, Vc, w, mu, k)
                
                return 
            end

            "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
            function thermalConductivityEstimate(
                Cp::Modelica.Media.Interfaces.Types.SpecificHeatCapacity,   # Constant pressure heat capacity
                eta::Modelica.Media.Interfaces.Types.DynamicViscosity,   # Dynamic viscosity
                method::Int64 = 1,               # 1: Eucken Method, 2: Modified Eucken Method
                data::IdealGases.Common.DataRecord,   # Ideal gas data
                # returns lambda::Modelica.Media.Interfaces.Types.ThermalConductivity,   # Thermal conductivity [W/(m.k)]
                )
                lambda = Functions.thermalConductivityEstimate(Cp, eta, method, data)
                
                return lambda
            end
        end
    end
end
end
