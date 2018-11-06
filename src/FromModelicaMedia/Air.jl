
"Medium models for air"
module Air
    #=
    @extends Modelica.Icons.VariantsPackage()
    =#

    "Air: Simple dry air model (0..100 degC)"
    module SimpleAir
        #=
        @extends Modelica.Icons.MaterialProperty()
        =#
        #=
        @extends Interfaces.PartialSimpleIdealGasMedium(
            mediumName="SimpleAir", 
            cp_const=1005.45, 
            MM_const=0.0289651159, 
            R_gas=Constants.R/0.0289651159, 
            eta_const=1.82e-5, 
            lambda_const=0.026, 
            T_min=Cv.from_degC(0), 
            T_max=Cv.from_degC(100), 
            fluidConstants=airConstants, 
            Temperature(
            min=Modelica.SIunits.Conversions.from_degC(0), 
            max=Modelica.SIunits.Conversions.from_degC(100)))
        =#
        const airConstants = Array{Modelica.Media.Interfaces.Types.Basic.FluidConstants, 1}=[Modelica.Media.Interfaces.Types.Basic.FluidConstants(iupacName="simple air", casRegistryNumber="not a real substance", chemicalFormula="N2, O2", structureFormula="N2, O2", molarMass=Modelica.Media.IdealGases.Common.SingleGasesData.N2.MM)]   # Constant data for the fluid
    end

    "Air: Detailed dry air model as ideal gas (200..6000 K)"
    module DryAirNasa
        #=
        @extends Modelica.Icons.MaterialProperty()
        =#
        #=
        @extends IdealGases.Common.SingleGasNasa(
            mediumName="Air", 
            data=IdealGases.Common.SingleGasesData.Air, 
            fluidConstants=[IdealGases.Common.FluidData.N2])
        =#

        "Return dynamic viscosity of dry air (simple polynomial, moisture influence small, valid from 123.15 K to 1273.15 K, outside of this range linear extrapolation is used)"
        function dynamicViscosity(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns eta::DynamicViscosity,   # Dynamic viscosity
            )
            eta = 1e-6*Polynomials_Temp.evaluateWithRange([9.7391102886305869E-15,-3.1353724870333906E-11,4.3004876595642225E-08,-3.8228016291758240E-05,5.0427874367180762E-02,1.7239260139242528E+01], Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T))
            
            return eta, 
        end

        "Return thermal conductivity of dry air (simple polynomial, moisture influence small, valid from 123.15 K to 1273.15 K, outside of this range linear extrapolation is used)"
        function thermalConductivity(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            method::Int64,                   # Dummy for compatibility reasons
            # returns lambda::ThermalConductivity,   # Thermal conductivity
            )
            lambda = 1e-3*Polynomials_Temp.evaluateWithRange([6.5691470817717812E-15,-3.4025961923050509E-11,5.3279284846303157E-08,-4.5340839289219472E-05,7.6129675309037664E-02,2.4169481088097051E+01], Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T))
            
            return lambda, , 
        end
    end

    "ReferenceAir: Detailed dry air model with a large operating range (130 ... 2000 K, 0 ... 2000 MPa) based on Helmholtz equations of state"
    module ReferenceAir
        #=
        @extends Modelica.Icons.VariantsPackage()
        =#
        const airConstants = Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants(
            chemicalFormula="N2+O2+Ar", 
            structureFormula="N2+O2+Ar", 
            casRegistryNumber="1", 
            iupacName="air", 
            molarMass=0.02896546, 
            criticalTemperature=132.5306, 
            criticalPressure=3.786e6, 
            criticalMolarVolume=0.02896546/342.68, 
            triplePointTemperature=63.05"From N2", 
            triplePointPressure=0.1253e5"From N2", 
            normalBoilingPoint=78.903, 
            meltingPoint=0, 
            acentricFactor=0.0335, 
            dipoleMoment=0.0, 
            hasCriticalData=true, 
            hasFundamentalEquation=true, 
            hasAccurateViscosityData=true, 
            hasAcentricFactor=true)

        "ReferenceAir.Air_ph: Detailed dry air model (130 ... 2000 K) explicit in p and h"
        module Air_ph
            #=
            @extends Modelica.Icons.MaterialProperty()
            =#
            #=
            @extends Modelica.Media.Air.ReferenceAir.Air_Base(
                ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph, 
                ph_explicit=true, 
                dT_explicit=false, 
                pT_explicit=false)
            =#
        end

        "ReferenceAir.Air_pT: Detailed dry air model (130 ... 2000 K) explicit in p and T"
        module Air_pT
            #=
            @extends Modelica.Icons.MaterialProperty()
            =#
            #=
            @extends Modelica.Media.Air.ReferenceAir.Air_Base(
                ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT, 
                ph_explicit=false, 
                dT_explicit=false, 
                pT_explicit=true)
            =#
        end

        "ReferenceAir.Air_dT: Detailed dry air model (130 ... 2000 K) explicit in d and T"
        module Air_dT
            #=
            @extends Modelica.Icons.MaterialProperty()
            =#
            #=
            @extends Modelica.Media.Air.ReferenceAir.Air_Base(
                ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.dTX, 
                ph_explicit=false, 
                dT_explicit=true, 
                pT_explicit=false)
            =#
        end

        "Properties of dry air calculated using the equation of state by Lemmon et. al."
        module Air_Base
            #=
            @extends Modelica.Media.Interfaces.PartialPureSubstance(
                mediumName="Air", 
                substanceNames=["air"], 
                singleState=false, 
                SpecificEnthalpy(
                start=1.0e5, 
                nominal=5.0e5), 
                Density(
                start=1.0, 
                nominal=1.2), 
                AbsolutePressure(
                start=1e5, 
                nominal=1e5, 
                min=1.0, 
                max=2000e6), 
                Temperature(
                start=273.15, 
                nominal=293.15, 
                min=130, 
                max=2000))
            =#
            const ph_explicit = Missing      # True if explicit in pressure and specific enthalpy
            const dT_explicit = Missing      # True if explicit in density and temperature
            const pT_explicit = Missing      # True if explicit in pressure and temperature
            mutable struct ThermodynamicState
                #= @extends ThermodynamicState
                 =#
                h::SpecificEnthalpy              # Specific enthalpy
                d::Density                       # Density
                T::Temperature                   # Temperature
                p::AbsolutePressure              # Pressure
            end
            @model BaseProperties begin
                #= @extends BaseProperties
                (
                h(
                stateSelect=if ph_explicit && preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
                d(
                stateSelect=if dT_explicit && preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
                T(
                stateSelect=if (pT_explicit || dT_explicit) && preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
                p(
                stateSelect=if (pT_explicit || ph_explicit) && preferredMediumStates; StateSelect.prefer else StateSelect.default end)) =#
                @inherits MM, ReferenceAir, dT_explicit, p, pressure_dT, d, T, h, specificEnthalpy_dT, ph_explicit, density_ph, temperature_ph, specificEnthalpy_pT, density_pT, u, R, Constants, state
                @equations begin
                                MM = ReferenceAir.Air_Utilities.Basic.Constants.MM
                                if dT_explicit
                    p = pressure_dT(d, T)
                    h = specificEnthalpy_dT(d, T)
                elseif ph_explicit
                    d = density_ph(p, h)
                    T = temperature_ph(p, h)
                else
                    h = specificEnthalpy_pT(p, T)
                    d = density_pT(p, T)
                end
                                u = h-p/d
                                R = Constants.R
                                h = state.h
                                p = state.p
                                T = state.T
                                d = state.d
                    end
                
            end

            "Computes density as a function of pressure and specific enthalpy"
            function density_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                h::SpecificEnthalpy,             # Specific enthalpy
                # returns d::Density,            # Density
                )
                d = Air_Utilities.rho_ph(p, h)
                
                return d
            end

            "Computes temperature as a function of pressure and specific enthalpy"
            function temperature_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                h::SpecificEnthalpy,             # Specific enthalpy
                # returns T::Temperature,        # Temperature
                )
                T = Air_Utilities.T_ph(p, h)
                
                return T
            end

            "Compute temperature from pressure and specific enthalpy"
            function temperature_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                s::SpecificEntropy,              # Specific entropy
                # returns T::Temperature,        # Temperature
                )
                T = Air_Utilities.T_ps(p, s)
                
                return T
            end

            "Computes density as a function of pressure and specific enthalpy"
            function density_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                s::SpecificEntropy,              # Specific entropy
                # returns d::Density,            # Density
                )
                d = Air_Utilities.rho_ps(p, s)
                
                return d
            end

            "Computes pressure as a function of density and temperature"
            function pressure_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Density,                      # Density
                T::Temperature,                  # Temperature
                # returns p::AbsolutePressure,   # Pressure
                )
                p = Air_Utilities.p_dT(d, T)
                
                return p
            end

            "Computes specific enthalpy as a function of density and temperature"
            function specificEnthalpy_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Density,                      # Density
                T::Temperature,                  # Temperature
                # returns h::SpecificEnthalpy,   # Specific enthalpy
                )
                h = Air_Utilities.h_dT(d, T)
                
                return h
            end

            "Computes specific enthalpy as a function of pressure and temperature"
            function specificEnthalpy_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                T::Temperature,                  # Temperature
                # returns h::SpecificEnthalpy,   # Specific enthalpy
                )
                h = Air_Utilities.h_pT(p, T)
                
                return h
            end

            "Computes specific enthalpy as a function of pressure and temperature"
            function specificEnthalpy_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                s::SpecificEntropy,              # Specific entropy
                # returns h::SpecificEnthalpy,   # Specific enthalpy
                )
                h = Air_Utilities.h_ps(p, s)
                
                return h
            end

            "Computes density as a function of pressure and temperature"
            function density_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::AbsolutePressure,             # Pressure
                T::Temperature,                  # Temperature
                # returns d::Density,            # Density
                )
                d = Air_Utilities.rho_pT(p, T)
                
                return d
            end
            function dynamicViscosity(
                #= @extends dynamicViscosity
                 =#
                )
                eta = Air_Utilities.Transport.eta_dT(state.d, state.T)
                
                return 
            end
            function thermalConductivity(
                #= @extends thermalConductivity
                 =#
                )
                lambda = Air_Utilities.Transport.lambda_dT(state.d, state.T)
                
                return 
            end
            function pressure(
                #= @extends pressure
                 =#
                )
                p = state.p
                
                return 
            end
            function temperature(
                #= @extends temperature
                 =#
                )
                T = state.T
                
                return 
            end
            function density(
                #= @extends density
                 =#
                )
                d = state.d
                
                return 
            end
            function specificEnthalpy(
                #= @extends specificEnthalpy
                 =#
                )
                h = state.h
                
                return 
            end
            function specificInternalEnergy(
                #= @extends specificInternalEnergy
                 =#
                )
                u = state.h-state.p/state.d
                
                return 
            end
            function specificGibbsEnergy(
                #= @extends specificGibbsEnergy
                 =#
                )
                g = state.h-state.T*specificEntropy(state)
                
                return 
            end
            function specificHelmholtzEnergy(
                #= @extends specificHelmholtzEnergy
                 =#
                )
                f = state.h-state.p/state.d-state.T*specificEntropy(state)
                
                return 
            end
            function specificEntropy(
                #= @extends specificEntropy
                 =#
                )
                if dT_explicit
                    s = Air_Utilities.s_dT(state.d, state.T)
                elseif pT_explicit
                    s = Air_Utilities.s_pT(state.p, state.T)
                else
                    s = Air_Utilities.s_ph(state.p, state.h)
                end
                
                return 
            end
            function specificHeatCapacityCp(
                #= @extends specificHeatCapacityCp
                 =#
                )
                if dT_explicit
                    cp = Air_Utilities.cp_dT(state.d, state.T)
                elseif pT_explicit
                    cp = Air_Utilities.cp_pT(state.p, state.T)
                else
                    cp = Air_Utilities.cp_ph(state.p, state.h)
                end
                
                return 
            end
            function specificHeatCapacityCv(
                #= @extends specificHeatCapacityCv
                 =#
                )
                if dT_explicit
                    cv = Air_Utilities.cv_dT(state.d, state.T)
                elseif pT_explicit
                    cv = Air_Utilities.cv_pT(state.p, state.T)
                else
                    cv = Air_Utilities.cv_ph(state.p, state.h)
                end
                
                return 
            end
            function isentropicExponent(
                #= @extends isentropicExponent
                 =#
                )
                if dT_explicit
                    gamma = Air_Utilities.isentropicExponent_dT(state.d, state.T)
                elseif pT_explicit
                    gamma = Air_Utilities.isentropicExponent_pT(state.p, state.T)
                else
                    gamma = Air_Utilities.isentropicExponent_ph(state.p, state.h)
                end
                
                return 
            end
            function isothermalCompressibility(
                #= @extends isothermalCompressibility
                 =#
                )
                if dT_explicit
                    kappa = Air_Utilities.kappa_dT(state.d, state.T)
                elseif pT_explicit
                    kappa = Air_Utilities.kappa_pT(state.p, state.T)
                else
                    kappa = Air_Utilities.kappa_ph(state.p, state.h)
                end
                
                return 
            end
            function isobaricExpansionCoefficient(
                #= @extends isobaricExpansionCoefficient
                 =#
                )
                if dT_explicit
                    beta = Air_Utilities.beta_dT(state.d, state.T)
                elseif pT_explicit
                    beta = Air_Utilities.beta_pT(state.p, state.T)
                else
                    beta = Air_Utilities.beta_ph(state.p, state.h)
                end
                
                return 
            end
            function velocityOfSound(
                #= @extends velocityOfSound
                 =#
                )
                if dT_explicit
                    a = Air_Utilities.velocityOfSound_dT(state.d, state.T)
                elseif pT_explicit
                    a = Air_Utilities.velocityOfSound_pT(state.p, state.T)
                else
                    a = Air_Utilities.velocityOfSound_ph(state.p, state.h)
                end
                
                return 
            end
            function density_derh_p(
                #= @extends density_derh_p
                 =#
                )
                ddhp = Air_Utilities.ddhp(state.p, state.h)
                
                return 
            end
            function density_derp_h(
                #= @extends density_derp_h
                 =#
                )
                ddph = Air_Utilities.ddph(state.p, state.h)
                
                return 
            end
            function setState_dTX(
                #= @extends setState_dTX
                 =#
                )
                state = ThermodynamicState(d=d, T=T, h=specificEnthalpy_dT(d, T), p=pressure_dT(d, T))
                
                return 
            end
            function setState_phX(
                #= @extends setState_phX
                 =#
                )
                state = ThermodynamicState(d=density_ph(p, h), T=temperature_ph(p, h), h=h, p=p)
                
                return 
            end
            function setState_psX(
                #= @extends setState_psX
                 =#
                )
                state = ThermodynamicState(d=density_ps(p, s), T=temperature_ps(p, s), h=specificEnthalpy_ps(p, s), p=p)
                
                return 
            end
            function setState_pTX(
                #= @extends setState_pTX
                 =#
                )
                state = ThermodynamicState(d=density_pT(p, T), T=T, h=specificEnthalpy_pT(p, T), p=p)
                
                return 
            end
            function setSmoothState(
                #= @extends setSmoothState
                 =#
                )
                state = ThermodynamicState(p=smoothStep(x, state_a.p, state_b.p, x_small), h=smoothStep(x, state_a.h, state_b.h, x_small), d=density_ph(smoothStep(x, state_a.p, state_b.p, x_small), smoothStep(x, state_a.h, state_b.h, x_small)), T=temperature_ph(smoothStep(x, state_a.p, state_b.p, x_small), smoothStep(x, state_a.h, state_b.h, x_small)))
                
                return 
            end
            function isentropicEnthalpy(
                #= @extends isentropicEnthalpy
                 =#
                )
                h_is = specificEnthalpy_psX(p_downstream, specificEntropy(refState), reference_X)
                
                return 
            end
            function molarMass(
                #= @extends molarMass
                 =#
                )
                MM = Modelica.Media.Air.ReferenceAir.airConstants.molarMass
                
                return 
            end
        end

        "Low level and utility computation for high accuracy dry air properties"
        module Air_Utilities
            #=
            @extends Modelica.Icons.UtilitiesPackage()
            =#

            "Fundamental equation of state"
            module Basic
                #=
                @extends Modelica.Icons.BasesPackage()
                =#
                const Constants = Modelica.Media.Common.FundamentalConstants(
                    R_bar=8.31451, 
                    R=287.117, 
                    MM=28.9586E-003, 
                    rhored=10447.7, 
                    Tred=132.6312, 
                    pred=3785020, 
                    h_off=1589557.62320524, 
                    s_off=6610.41237132543)

                "Helmholtz equation of state"
                function Helmholtz(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature (K)
                    # returns f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and derivatives w.r.t. delta and tau
                    )
                    const N_0::Array{Float64, 1} = Missing,
                    const N::Array{Float64, 1} = Missing,
                    const i::Array{Int64, 1} = Missing,
                    const j::Array{Float64, 1} = Missing,
                    const l::Array{Int64, 1} = Missing,
                    f.d = d
                                        f.T = T
                                        f.R = ReferenceAir.Air_Utilities.Basic.Constants.R
                                        f.delta = d/
                        (ReferenceAir.Air_Utilities.Basic.Constants.MM*ReferenceAir.Air_Utilities.Basic.Constants.rhored)
                                        f.tau = ReferenceAir.Air_Utilities.Basic.Constants.Tred/T
                                        f.f = 0
                                        for k in 1:5
                    f.f = f.f+N_0[k]*f.tau^(k-4)
                    end
                    
                                        f.f = f.f+log(f.delta)+N_0[6]*f.tau*sqrt(f.tau)+N_0[7]*log(f.tau)+N_0[8]*log(1-exp(-N_0[11]*f.tau))+N_0[9]*log(1-exp(-N_0[12]*f.tau))+N_0[10]*log(2/3+exp(N_0[13]*f.tau))
                                        for k in 1:10
                    f.f = f.f+N[k]*f.delta^i[k]*f.tau^j[k]
                    end
                    
                                        for k in 11:19
                    f.f = f.f+N[k]*f.delta^i[k]*f.tau^j[k]*exp(-f.delta^l[k])
                    end
                    
                                        f.fdelta = 0
                                        f.fdelta = 1/f.delta
                                        for k in 1:10
                    f.fdelta = f.fdelta+i[k]*N[k]*f.delta^(i[k]-1)*f.tau^j[k]
                    end
                    
                                        for k in 11:19
                    f.fdelta = f.fdelta+N[k]*f.delta^(i[k]-1)*f.tau^j[k]*exp(-f.delta^l[k])*(i[k]-l[k]*f.delta^l[k])
                    end
                    
                                        f.fdeltadelta = 0
                                        f.fdeltadelta = -1/f.delta^2
                                        for k in 1:10
                    f.fdeltadelta = f.fdeltadelta+i[k]*(i[k]-1)*N[k]*f.delta^(i[k]-2)*f.tau^j[k]
                    end
                    
                                        for k in 11:19
                    f.fdeltadelta = f.fdeltadelta+N[k]*f.delta^(i[k]-2)*f.tau^j[k]*exp(-f.delta^l[k])*
                        ((i[k]-l[k]*f.delta^l[k])*(i[k]-1-l[k]*f.delta^l[k])-l[k]^2*f.delta^l[k])
                    end
                    
                                        f.ftau = 0
                                        for k in 1:5
                    f.ftau = f.ftau+(k-4)*N_0[k]*f.tau^(k-5)
                    end
                    
                                        f.ftau = f.ftau+1.5*N_0[6]*sqrt(f.tau)+N_0[7]/f.tau+N_0[8]*N_0[11]/(exp(N_0[11]*f.tau)-1)+N_0[9]*N_0[12]/(exp(N_0[12]*f.tau)-1)+N_0[10]*N_0[13]/(2/3*exp(-N_0[13]*f.tau)+1)
                                        for k in 1:10
                    f.ftau = f.ftau+j[k]*N[k]*f.delta^i[k]*f.tau^(j[k]-1)
                    end
                    
                                        for k in 11:19
                    f.ftau = f.ftau+j[k]*N[k]*f.delta^i[k]*f.tau^(j[k]-1)*exp(-f.delta^l[k])
                    end
                    
                                        f.ftautau = 0
                                        for k in 1:3
                    f.ftautau = f.ftautau+(k-4)*(k-5)*N_0[k]*f.tau^(k-6)
                    end
                    
                                        f.ftautau = f.ftautau+0.75*N_0[6]/sqrt(f.tau)-N_0[7]/f.tau^2-N_0[8]*N_0[11]^2*exp(N_0[11]*f.tau)/(exp(N_0[11]*f.tau)-1)^2-N_0[9]*N_0[12]^2*exp(N_0[12]*f.tau)/(exp(N_0[12]*f.tau)-1)^2+2/3*N_0[10]*N_0[13]^2*exp(-N_0[13]*f.tau)/(2/3*exp(-N_0[13]*f.tau)+1)^2
                                        for k in 1:10
                    f.ftautau = f.ftautau+j[k]*(j[k]-1)*N[k]*f.delta^i[k]*f.tau^(j[k]-2)
                    end
                    
                                        for k in 11:19
                    f.ftautau = f.ftautau+j[k]*(j[k]-1)*N[k]*f.delta^i[k]*f.tau^(j[k]-2)*exp(-f.delta^l[k])
                    end
                    
                                        f.fdeltatau = 0
                                        for k in 1:10
                    f.fdeltatau = f.fdeltatau+i[k]*j[k]*N[k]*f.delta^(i[k]-1)*f.tau^(j[k]-1)
                    end
                    
                                        for k in 11:19
                    f.fdeltatau = f.fdeltatau+j[k]*N[k]*f.delta^(i[k]-1)*f.tau^(j[k]-1)*exp(-f.delta^l[k])*(i[k]-l[k]*f.delta^l[k])
                    end
                    
                    
                    return f
                end
            end

            "Inverse function"
            module Inverses
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Accuracy of the iterations"
                mutable struct accuracy
                    #=
                    @extends Modelica.Icons.Record()
                    =#
                    const delp::Float64 = Missing    # Accuracy of p
                    const delh::Float64 = Missing    # Accuracy of h
                    const dels::Float64 = Missing    # Accuracy of s
                end

                "Compute d for given p and T"
                function dofpT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature (K)
                    delp::Float64,                   # Iteration converged if (p-pre(p) < delp)
                    # returns d::Float64,            # Density
                    )
                    i::Int64,                        # Loop counter
                    dp::Float64,                     # Pressure difference
                    deld::Float64,                   # Density step
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    nDerivs::Modelica.Media.Common.NewtonDerivatives_pT,   # Derivatives needed in Newton iteration
                    found::Bool,                     # Flag for iteration success
                    d = p/
                        (ReferenceAir.Air_Utilities.Basic.Constants.R*T)
                                        while ((i<100) && ! found)
                    f = Basic.Helmholtz(d, T)
                    nDerivs = Modelica.Media.Common.Helmholtz_pT(f)
                    dp = nDerivs.p-p
                    if (abs(dp)<=delp)
                        found = true
                    end
                    deld = dp/nDerivs.pd
                    d = d-deld
                    i = i+1
                    end
                    
                    
                    return d
                end

                "Return d and T as a function of p and h"
                function dTofph(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    h::Float64,                      # Specific enthalpy
                    delp::Float64,                   # Iteration accuracy
                    delh::Float64,                   # Iteration accuracy
                    # returns d::Float64,            # Density
                    # returns T::Float64,            # Temperature (K)
                    )
                    Tguess::Float64,                 # Initial temperature
                    dguess::Float64,                 # Initial density
                    i::Int64,                        # Iteration counter
                    dh::Float64,                     # Newton-error in h-direction
                    dp::Float64,                     # Newton-error in p-direction
                    det::Float64,                    # Determinant of directional derivatives
                    deld::Float64,                   # Newton-step in d-direction
                    delt::Float64,                   # Newton-step in T-direction
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    nDerivs::Modelica.Media.Common.NewtonDerivatives_ph,   # Derivatives needed in Newton iteration
                    found::Bool,                     # Flag for iteration success
                    T = h/1000+273.15
                                        d = p/
                        (ReferenceAir.Air_Utilities.Basic.Constants.R*T)
                                        i = 0
                                        while ((i<100) && ! found)
                    f = Basic.Helmholtz(d, T)
                    nDerivs = Modelica.Media.Common.Helmholtz_ph(f)
                    dh = nDerivs.h-ReferenceAir.Air_Utilities.Basic.Constants.h_off-h
                    dp = nDerivs.p-p
                    if ((abs(dh)<=delh) && (abs(dp)<=delp))
                        found = true
                    end
                    det = nDerivs.ht*nDerivs.pd-nDerivs.pt*nDerivs.hd
                    delt = (nDerivs.pd*dh-nDerivs.hd*dp)/det
                    deld = (nDerivs.ht*dp-nDerivs.pt*dh)/det
                    T = T-delt
                    d = d-deld
                    i = i+1
                    end
                    
                    
                    return d, T
                end

                "Return d and T as a function of p and s"
                function dTofps(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    s::Float64,                      # Specific entropy
                    delp::Float64,                   # Iteration accuracy
                    dels::Float64,                   # Iteration accuracy
                    # returns d::Float64,            # Density
                    # returns T::Float64,            # Temperature (K)
                    )
                    Tguess::Float64,                 # Initial temperature
                    dguess::Float64,                 # Initial density
                    i::Int64,                        # Iteration counter
                    ds::Float64,                     # Newton-error in s-direction
                    dp::Float64,                     # Newton-error in p-direction
                    det::Float64,                    # Determinant of directional derivatives
                    deld::Float64,                   # Newton-step in d-direction
                    delt::Float64,                   # Newton-step in T-direction
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    nDerivs::Modelica.Media.Common.NewtonDerivatives_ps,   # Derivatives needed in Newton iteration
                    found::Bool,                     # Flag for iteration success
                    T = 273.15
                                        d = p/
                        (ReferenceAir.Air_Utilities.Basic.Constants.R*T)
                                        i = 0
                                        while ((i<100) && ! found)
                    f = Basic.Helmholtz(d, T)
                    nDerivs = Modelica.Media.Common.Helmholtz_ps(f)
                    ds = nDerivs.s-ReferenceAir.Air_Utilities.Basic.Constants.s_off-s
                    dp = nDerivs.p-p
                    if ((abs(ds)<=dels) && (abs(dp)<=delp))
                        found = true
                    end
                    det = nDerivs.st*nDerivs.pd-nDerivs.pt*nDerivs.sd
                    delt = (nDerivs.pd*ds-nDerivs.sd*dp)/det
                    deld = (nDerivs.st*dp-nDerivs.pt*ds)/det
                    T = T-delt
                    d = d-deld
                    i = i+1
                    end
                    
                    
                    return d, T
                end
            end

            "Transport properties for air"
            module Transport
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Return dynamic viscosity as a function of d and T"
                function eta_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns eta::Float64,          # Dynamic viscosity
                    )
                    delta::Float64,                  # Reduced density
                    tau::Float64,                    # Reciprocal reduced temperature
                    Omega::Float64,                  # Collision integral
                    eta_0::Float64,                  # Dilute gas viscosity
                    eta_r::Float64,                  # Residual fluid viscosity
                    const b::Array{Float64, 1} = Missing,
                    const Nvis::Array{Float64, 1} = Missing,
                    const tvis::Array{Float64, 1} = Missing,
                    const dvis::Array{Int64, 1} = Missing,
                    const lvis::Array{Int64, 1} = Missing,
                    const gammavis::Array{Int64, 1} = Missing,
                    Omega = exp(Modelica.Media.Incompressible.TableBased.Polynomials_Temp.evaluate([b[5],b[4],b[3],b[2],b[1]], log(T/103.3)))
                                        eta_0 = 0.0266958*sqrt(1000*ReferenceAir.Air_Utilities.Basic.Constants.MM*T)/(0.36^2*Omega)
                                        for i in 1:5
                    eta_r = eta_r+
                        (Nvis[i]*(tau^tvis[i])*(delta^dvis[i])*exp(-gammavis[i]*(delta^lvis[i])))
                    end
                    
                                        eta = (eta_0+eta_r)*1E-006
                    
                    return eta
                end

                "Return thermal conductivity as a function of d and T"
                function lambda_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns lambda::Float64,       # Thermal conductivity
                    )
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    lambda_0::Float64,               # Dilute gas thermal conductivity
                    lambda_r::Float64,               # Residual fluid thermal conductivity
                    lambda_c::Float64,               # Thermal conductivity critical enhancement
                    Omega::Float64,                  # Collision integral
                    eta_0::Float64,                  # Dilute gas viscosity
                    pddT::Float64,                
                    pddTref::Float64,             
                    pdTp::Float64,                
                    xi::Float64,                  
                    xiref::Float64,               
                    Omega_tilde::Float64,         
                    Omega_0_tilde::Float64,       
                    cv::Float64,                  
                    cp::Float64,                  
                    const b::Array{Float64, 1} = Missing,
                    const Ncon::Array{Float64, 1} = Missing,
                    const tcon::Array{Float64, 1} = Missing,
                    const dcon::Array{Int64, 1} = Missing,
                    const lcon::Array{Int64, 1} = Missing,
                    const gammacon::Array{Int64, 1} = Missing,
                    f = Basic.Helmholtz(d, 265.262)
                                        pddTref = ReferenceAir.Air_Utilities.Basic.Constants.R_bar*265.262*
                        (1+2*f.delta*(f.fdelta-1/f.delta)+f.delta^2*(f.fdeltadelta+1/f.delta^2))
                                        xiref = ReferenceAir.Air_Utilities.Basic.Constants.pred*
                        (d/ReferenceAir.Air_Utilities.Basic.Constants.MM)/ReferenceAir.Air_Utilities.Basic.Constants.rhored^2/pddTref
                                        f = Basic.Helmholtz(d, T)
                                        Omega = exp(Modelica.Media.Incompressible.TableBased.Polynomials_Temp.evaluate([b[5],b[4],b[3],b[2],b[1]], log(T/103.3)))
                                        eta_0 = 0.0266958*sqrt(1000*ReferenceAir.Air_Utilities.Basic.Constants.MM*T)/(0.36^2*Omega)
                                        lambda_0 = Ncon[1]*eta_0+Ncon[2]*f.tau^tcon[2]+Ncon[3]*f.tau^tcon[3]
                                        for i in 4:9
                    lambda_r = lambda_r+Ncon[i]*f.tau^tcon[i]*f.delta^dcon[i]*exp(-gammacon[i]*f.delta^lcon[i])
                    end
                    
                                        pddT = ReferenceAir.Air_Utilities.Basic.Constants.R*T*
                        (1+2*f.delta*(f.fdelta-1/f.delta)+f.delta^2*(f.fdeltadelta+1/f.delta^2))
                                        xi = ReferenceAir.Air_Utilities.Basic.Constants.pred*
                        (d/ReferenceAir.Air_Utilities.Basic.Constants.MM)/ReferenceAir.Air_Utilities.Basic.Constants.rhored^2/
                        (pddT*ReferenceAir.Air_Utilities.Basic.Constants.MM)
                                        xi = xi-xiref*265.262/T
                                        if (xi<=0)
                        lambda_c = 0
                    else
                        xi = 0.11*(xi/0.055)^(0.63/1.2415)
                        pdTp = ReferenceAir.Air_Utilities.Basic.Constants.R*d*
                        (1+f.delta*(f.fdelta-1/f.delta)-f.delta*f.tau*f.fdeltatau)
                        cv = ReferenceAir.Air_Utilities.Basic.Constants.R*(-f.tau*f.tau*f.ftautau)
                        cp = cv+T*pdTp*pdTp/(d*d*pddT)
                        Omega_tilde = 2/Modelica.Constants.pi*((cp-cv)/cp*atan(xi/0.31)+cv/cp*xi/0.31)
                        Omega_0_tilde = 2/Modelica.Constants.pi*
                        (1-exp(-1/
                        ((0.31/xi)+1/3*(xi/0.31)^2*
                        (ReferenceAir.Air_Utilities.Basic.Constants.rhored/
                        (d/ReferenceAir.Air_Utilities.Basic.Constants.MM))^2)))
                        lambda_c = d*cp*1.380658E-023*1.01*T/
                        (6*Modelica.Constants.pi*xi*eta_dT(d, T))*(Omega_tilde-Omega_0_tilde)*1E012
                    end
                                        lambda = (lambda_0+lambda_r+lambda_c)/1000
                    
                    return lambda
                end
            end

            "Intermediate property record for air"
            function airBaseProp_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                # returns aux::Common.AuxiliaryProperties,   # Auxiliary record
                )
                f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                aux.p = p
                                aux.s = s
                                aux.R = ReferenceAir.Air_Utilities.Basic.Constants.R
                                (aux.rhoaux.T) = Inverses.dTofps(p=p, s=s, delp=iter.delp, dels=iter.dels)
                                f = Basic.Helmholtz(aux.rho, aux.T)
                                aux.h = aux.R*aux.T*(f.tau*f.ftau+f.delta*f.fdelta)-ReferenceAir.Air_Utilities.Basic.Constants.h_off
                                aux.pd = aux.R*aux.T*f.delta*(2*f.fdelta+f.delta*f.fdeltadelta)
                                aux.pt = aux.R*aux.rho*f.delta*(f.fdelta-f.tau*f.fdeltatau)
                                aux.cv = aux.R*(-f.tau*f.tau*f.ftautau)
                                aux.cp = aux.cv+aux.T*aux.pt*aux.pt/(aux.rho*aux.rho*aux.pd)
                                aux.vp = -1/(aux.rho*aux.rho)*1/aux.pd
                                aux.vt = aux.pt/(aux.rho*aux.rho*aux.pd)
                
                return aux
            end

            "Density as function of pressure and specific entropy"
            function rho_props_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns rho::Float64,          # Density
                )
                rho = aux.rho
                
                return rho
            end

            "Density as function of pressure and specific entropy"
            function rho_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                # returns rho::Float64,          # Density
                )
                rho = rho_props_ps(p, s, Air_Utilities.airBaseProp_ps(p, s))
                
                return rho
            end

            "Temperature as function of pressure and specific entropy"
            function T_props_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns T::Float64,            # Temperature
                )
                T = aux.T
                
                return T
            end

            "Temperature as function of pressure and specific entropy"
            function T_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                # returns T::Float64,            # Temperature
                )
                T = T_props_ps(p, s, Air_Utilities.airBaseProp_ps(p, s))
                
                return T
            end

            "Specific enthalpy as function or pressure and temperature"
            function h_props_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns h::Float64,            # Specific enthalpy
                )
                h = aux.h
                
                return h
            end

            "Specific enthalpy as function or pressure and temperature"
            function h_ps(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                s::Float64,                      # Specific entropy
                # returns h::Float64,            # Specific enthalpy
                )
                h = h_props_ps(p, s, Air_Utilities.airBaseProp_ps(p, s))
                
                return h
            end

            "Intermediate property record for air"
            function airBaseProp_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns aux::Common.AuxiliaryProperties,   # Auxiliary record
                )
                f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                error::Int64,                    # Error flag for inverse iterations
                aux.p = p
                                aux.h = h
                                aux.R = ReferenceAir.Air_Utilities.Basic.Constants.R
                                (aux.rhoaux.T) = Inverses.dTofph(p, h, delp=iter.delp, delh=iter.delh)
                                f = Basic.Helmholtz(aux.rho, aux.T)
                                aux.s = aux.R*(f.tau*f.ftau-f.f)-ReferenceAir.Air_Utilities.Basic.Constants.s_off
                                aux.pd = aux.R*aux.T*f.delta*(2*f.fdelta+f.delta*f.fdeltadelta)
                                aux.pt = aux.R*aux.rho*f.delta*(f.fdelta-f.tau*f.fdeltatau)
                                aux.cv = aux.R*(-f.tau*f.tau*f.ftautau)
                                aux.cp = aux.cv+aux.T*aux.pt*aux.pt/(aux.rho*aux.rho*aux.pd)
                                aux.vp = -1/(aux.rho*aux.rho)*1/aux.pd
                                aux.vt = aux.pt/(aux.rho*aux.rho*aux.pd)
                
                return aux
            end

            "Density as function of pressure and specific enthalpy"
            function rho_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns rho::Float64,          # Density
                )
                rho = aux.rho
                
                return rho
            end

            "Density as function of pressure and specific enthalpy"
            function rho_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns rho::Float64,          # Density
                )
                rho = rho_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return rho
            end

            "Derivative function of rho_ph"
            function rho_ph_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                p_der::Float64,                  # Derivative of pressure
                h_der::Float64,                  # Derivative of specific enthalpy
                # returns rho_der::Float64,      # Derivative of density
                )
                rho_der = 
                    ((aux.rho*(aux.cv*aux.rho+aux.pt))/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt))*p_der+
                    (-aux.rho*aux.rho*aux.pt/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt))*h_der
                
                return rho_der
            end

            "Temperature as function of pressure and specific enthalpy"
            function T_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns T::Float64,            # Temperature
                )
                T = aux.T
                
                return T
            end

            "Temperature as function of pressure and specific enthalpy"
            function T_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns T::Float64,            # Temperature
                )
                T = T_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return T
            end

            "Derivative function of T_ph"
            function T_ph_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                p_der::Float64,                  # Derivative of pressure
                h_der::Float64,                  # Derivative of specific enthalpy
                # returns T_der::Float64,        # Derivative of temperature
                )
                T_der = 
                    ((-aux.rho*aux.pd+aux.T*aux.pt)/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt))*p_der+
                    ((aux.rho*aux.rho*aux.pd)/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt))*h_der
                
                return T_der
            end

            "Specific entropy as function of pressure and specific enthalpy"
            function s_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns s::Float64,            # Specific entropy
                )
                s = aux.s
                
                return s
            end

            "Specific entropy as function of pressure and specific enthalpy"
            function s_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns s::Float64,            # Specific entropy
                )
                s = s_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return s
            end

            "Specific entropy as function of pressure and specific enthalpy"
            function s_ph_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                p_der::Float64,                  # Derivative of pressure
                h_der::Float64,                  # Derivative of specific enthalpy
                # returns s_der::Float64,        # Derivative of entropy
                )
                s_der = -1/(aux.rho*aux.T)*p_der+1/aux.T*h_der
                
                return s_der
            end

            "Specific heat capacity at constant volume as function of pressure and specific enthalpy"
            function cv_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = aux.cv
                
                return cv
            end

            "Specific heat capacity at constant volume as function of pressure and specific enthalpy"
            function cv_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = cv_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return cv
            end

            "Specific heat capacity at constant pressure as function of pressure and specific enthalpy"
            function cp_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = aux.cp
                
                return cp
            end

            "Specific heat capacity at constant pressure as function of pressure and specific enthalpy"
            function cp_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = cp_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return cp
            end

            "Isobaric expansion coefficient as function of pressure and specific enthalpy"
            function beta_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = aux.pt/(aux.rho*aux.pd)
                
                return beta
            end

            "Isobaric expansion coefficient as function of pressure and specific enthalpy"
            function beta_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = beta_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return beta
            end

            "Isothermal compressibility factor as function of pressure and specific enthalpy"
            function kappa_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = 1/(aux.rho*aux.pd)
                
                return kappa
            end

            "Isothermal compressibility factor as function of pressure and specific enthalpy"
            function kappa_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = kappa_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return kappa
            end

            "Speed of sound as function of pressure and specific enthalpy"
            function velocityOfSound_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns a::Float64,            # Speed of sound
                )
                a = sqrt(max(0, aux.pd+aux.pt*aux.pt*aux.T/(aux.rho*aux.rho*aux.cv)))
                
                return a
            end

            
            function velocityOfSound_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns a::Float64,            # Speed of sound
                )
                a = velocityOfSound_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return a
            end

            "Isentropic exponent as function of pressure and specific enthalpy"
            function isentropicExponent_props_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = 1/(aux.rho*p)*
                    (
                    (aux.pd*aux.cv*aux.rho*aux.rho+aux.pt*aux.pt*aux.T)/(aux.cv))
                
                return gamma
            end

            "Isentropic exponent as function of pressure and specific enthalpy"
            function isentropicExponent_ph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = isentropicExponent_props_ph(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return gamma
            end

            "Density derivative by pressure"
            function ddph_props(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns ddph::Float64,         # Density derivative by pressure
                )
                ddph = 
                    ((aux.rho*(aux.cv*aux.rho+aux.pt))/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt))
                
                return ddph
            end

            "Density derivative by pressure"
            function ddph(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns ddph::Float64,         # Density derivative by pressure
                )
                ddph = ddph_props(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return ddph
            end

            "Density derivative by specific enthalpy"
            function ddhp_props(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns ddhp::Float64,         # Density derivative by specific enthalpy
                )
                ddhp = -aux.rho*aux.rho*aux.pt/
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt)
                
                return ddhp
            end

            "Density derivative by specific enthalpy"
            function ddhp(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                h::Float64,                      # Specific enthalpy
                # returns ddhp::Float64,         # Density derivative by specific enthalpy
                )
                ddhp = ddhp_props(p, h, Air_Utilities.airBaseProp_ph(p, h))
                
                return ddhp
            end

            "Intermediate property record for air (p and T preferred states)"
            function airBaseProp_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns aux::Common.AuxiliaryProperties,   # Auxiliary record
                )
                f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                aux.p = p
                                aux.T = T
                                aux.R = ReferenceAir.Air_Utilities.Basic.Constants.R
                                (aux.rho) = Inverses.dofpT(p=p, T=T, delp=iter.delp)
                                f = Basic.Helmholtz(aux.rho, T)
                                aux.h = aux.R*T*(f.tau*f.ftau+f.delta*f.fdelta)-ReferenceAir.Air_Utilities.Basic.Constants.h_off
                                aux.s = aux.R*(f.tau*f.ftau-f.f)-ReferenceAir.Air_Utilities.Basic.Constants.s_off
                                aux.pd = aux.R*T*f.delta*(2*f.fdelta+f.delta*f.fdeltadelta)
                                aux.pt = aux.R*aux.rho*f.delta*(f.fdelta-f.tau*f.fdeltatau)
                                aux.cv = aux.R*(-f.tau*f.tau*f.ftautau)
                                aux.cp = aux.cv+aux.T*aux.pt*aux.pt/(aux.rho*aux.rho*aux.pd)
                                aux.vp = -1/(aux.rho*aux.rho)*1/aux.pd
                                aux.vt = aux.pt/(aux.rho*aux.rho*aux.pd)
                
                return aux
            end

            "Density as function or pressure and temperature"
            function rho_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns rho::Float64,          # Density
                )
                rho = aux.rho
                
                return rho
            end

            "Density as function or pressure and temperature"
            function rho_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns rho::Float64,          # Density
                )
                rho = rho_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return rho
            end

            "Derivative function of rho_pT"
            function rho_pT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                # returns rho_der::Float64,      # Derivative of density
                )
                rho_der = (1/aux.pd)*p_der-(aux.pt/aux.pd)*T_der
                
                return rho_der
            end

            "Specific enthalpy as function or pressure and temperature"
            function h_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns h::Float64,            # Specific enthalpy
                )
                h = aux.h
                
                return h
            end

            "Specific enthalpy as function or pressure and temperature"
            function h_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns h::Float64,            # Specific enthalpy
                )
                h = h_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return h
            end

            "Derivative function of h_pT"
            function h_pT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                # returns h_der::Float64,        # Derivative of specific enthalpy
                )
                h_der = 
                    ((-aux.rho*aux.pd+T*aux.pt)/(aux.rho*aux.rho*aux.pd))*p_der+
                    (
                    (aux.rho*aux.rho*aux.pd*aux.cv+aux.T*aux.pt*aux.pt)/(aux.rho*aux.rho*aux.pd))*T_der
                
                return h_der
            end

            "Specific entropy as function of pressure and temperature"
            function s_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns s::Float64,            # Specific entropy
                )
                s = aux.s
                
                return s
            end

            "Temperature as function of pressure and temperature"
            function s_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns s::Float64,            # Specific entropy
                )
                s = s_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return s
            end

            "Specific heat capacity at constant volume as function of pressure and temperature"
            function cv_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = aux.cv
                
                return cv
            end

            "Specific heat capacity at constant volume as function of pressure and temperature"
            function cv_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = cv_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return cv
            end

            "Specific heat capacity at constant pressure as function of pressure and temperature"
            function cp_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = aux.cp
                
                return cp
            end

            "Specific heat capacity at constant pressure as function of pressure and temperature"
            function cp_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = cp_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return cp
            end

            "Isobaric expansion coefficient as function of pressure and temperature"
            function beta_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = aux.pt/(aux.rho*aux.pd)
                
                return beta
            end

            "Isobaric expansion coefficient as function of pressure and temperature"
            function beta_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = beta_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return beta
            end

            "Isothermal compressibility factor as function of pressure and temperature"
            function kappa_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = 1/(aux.rho*aux.pd)
                
                return kappa
            end

            "Isothermal compressibility factor as function of pressure and temperature"
            function kappa_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = kappa_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return kappa
            end

            "Speed of sound as function of pressure and temperature"
            function velocityOfSound_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns a::Float64,            # Speed of sound
                )
                a = sqrt(max(0, 
                    (aux.pd*aux.rho*aux.rho*aux.cv+aux.pt*aux.pt*aux.T)/(aux.rho*aux.rho*aux.cv)))
                
                return a
            end

            "Speed of sound as function of pressure and temperature"
            function velocityOfSound_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns a::Float64,            # Speed of sound
                )
                a = velocityOfSound_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return a
            end

            "Isentropic exponent as function of pressure and temperature"
            function isentropicExponent_props_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = 1/(aux.rho*p)*
                    (
                    (aux.pd*aux.cv*aux.rho*aux.rho+aux.pt*aux.pt*aux.T)/(aux.cv))
                
                return gamma
            end

            "Isentropic exponent as function of pressure and temperature"
            function isentropicExponent_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = isentropicExponent_props_pT(p, T, Air_Utilities.airBaseProp_pT(p, T))
                
                return gamma
            end

            "Intermediate property record for air (d and T preferred states)"
            function airBaseProp_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns aux::Common.AuxiliaryProperties,   # Auxiliary record
                )
                f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                aux.rho = d
                                aux.T = T
                                aux.R = ReferenceAir.Air_Utilities.Basic.Constants.R
                                f = Basic.Helmholtz(d, T)
                                aux.p = aux.R*d*T*f.delta*f.fdelta
                                aux.h = aux.R*T*(f.tau*f.ftau+f.delta*f.fdelta)-ReferenceAir.Air_Utilities.Basic.Constants.h_off
                                aux.s = aux.R*(f.tau*f.ftau-f.f)-ReferenceAir.Air_Utilities.Basic.Constants.s_off
                                aux.pd = aux.R*T*f.delta*(2*f.fdelta+f.delta*f.fdeltadelta)
                                aux.pt = aux.R*d*f.delta*(f.fdelta-f.tau*f.fdeltatau)
                                aux.cv = aux.R*(-f.tau*f.tau*f.ftautau)
                                aux.cp = aux.cv+aux.T*aux.pt*aux.pt/(d*d*aux.pd)
                                aux.vp = -1/(aux.rho*aux.rho)*1/aux.pd
                                aux.vt = aux.pt/(aux.rho*aux.rho*aux.pd)
                
                return aux
            end

            "Specific enthalpy as function of density and temperature"
            function h_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns h::Float64,            # Specific enthalpy
                )
                h = aux.h
                
                return h
            end

            "Specific enthalpy as function of density and temperature"
            function h_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns h::Float64,            # Specific enthalpy
                )
                h = h_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return h
            end

            "Derivative function of h_dT"
            function h_dT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                d_der::Float64,                  # Derivative of density
                T_der::Float64,                  # Derivative of temperature
                # returns h_der::Float64,        # Derivative of specific enthalpy
                )
                h_der = ((-d*aux.pd+T*aux.pt)/(d*d))*d_der+((aux.cv*d+aux.pt)/d)*T_der
                
                return h_der
            end

            "Pressure as function of density and temperature"
            function p_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns p::Float64,            # Pressure
                )
                p = aux.p
                
                return p
            end

            "Pressure as function of density and temperature"
            function p_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns p::Float64,            # Pressure
                )
                p = p_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return p
            end

            "Derivative function of p_dT"
            function p_dT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                d_der::Float64,                  # Derivative of density
                T_der::Float64,                  # Derivative of temperature
                # returns p_der::Float64,        # Derivative of pressure
                )
                p_der = aux.pd*d_der+aux.pt*T_der
                
                return p_der
            end

            "Specific entropy as function of density and temperature"
            function s_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns s::Float64,            # Specific entropy
                )
                s = aux.s
                
                return s
            end

            "Temperature as function of density and temperature"
            function s_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns s::Float64,            # Specific entropy
                )
                s = s_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return s
            end

            "Specific heat capacity at constant volume as function of density and temperature"
            function cv_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = aux.cv
                
                return cv
            end

            "Specific heat capacity at constant volume as function of density and temperature"
            function cv_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns cv::Float64,           # Specific heat capacity
                )
                cv = cv_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return cv
            end

            "Specific heat capacity at constant pressure as function of density and temperature"
            function cp_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = aux.cp
                
                return cp
            end

            "Specific heat capacity at constant pressure as function of density and temperature"
            function cp_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns cp::Float64,           # Specific heat capacity
                )
                cp = cp_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return cp
            end

            "Isobaric expansion coefficient as function of density and temperature"
            function beta_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = aux.pt/(aux.rho*aux.pd)
                
                return beta
            end

            "Isobaric expansion coefficient as function of density and temperature"
            function beta_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns beta::Float64,         # Isobaric expansion coefficient
                )
                beta = beta_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return beta
            end

            "Isothermal compressibility factor as function of density and temperature"
            function kappa_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = 1/(aux.rho*aux.pd)
                
                return kappa
            end

            "Isothermal compressibility factor as function of density and temperature"
            function kappa_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns kappa::Float64,        # Isothermal compressibility factor
                )
                kappa = kappa_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return kappa
            end

            "Speed of sound as function of density and temperature"
            function velocityOfSound_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns a::Float64,            # Speed of sound
                )
                a = sqrt(max(0, 
                    (
                    (aux.pd*aux.rho*aux.rho*aux.cv+aux.pt*aux.pt*aux.T)/(aux.rho*aux.rho*aux.cv))))
                
                return a
            end

            "Speed of sound as function of density and temperature"
            function velocityOfSound_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns a::Float64,            # Speed of sound
                )
                a = velocityOfSound_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return a
            end

            "Isentropic exponent as function of density and temperature"
            function isentropicExponent_props_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                aux::Common.AuxiliaryProperties,   # Auxiliary record
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = 1/(aux.rho*aux.p)*
                    (
                    (aux.pd*aux.cv*aux.rho*aux.rho+aux.pt*aux.pt*aux.T)/(aux.cv))
                
                return gamma
            end

            "Isentropic exponent as function of density and temperature"
            function isentropicExponent_dT(
                #=
                @extends Modelica.Icons.Function()
                =#
                d::Float64,                      # Density
                T::Float64,                      # Temperature
                # returns gamma::Float64,        # Isentropic exponent
                )
                gamma = isentropicExponent_props_dT(d, T, Air_Utilities.airBaseProp_dT(d, T))
                
                return gamma
            end

            
            module ThermoFluidSpecial

                "Calculate the property record for dynamic simulation properties using p,h as states"
                function air_ph(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    h::Float64,                      # Specific enthalpy
                    # returns pro::Modelica.Media.Common.ThermoFluidSpecial.ThermoProperties_ph,   # Property record for dynamic simulation
                    )
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    T::Float64,                      # Temperature
                    d::Float64,                      # Density
                    (dT) = Air_Utilities.Inverses.dTofph(p=p, h=h, delp=1.0e-7, delh=1.0e-6)
                                        f = Air_Utilities.Basic.Helmholtz(d, T)
                                        pro = Modelica.Media.Common.ThermoFluidSpecial.helmholtzToProps_ph(f)
                    
                    return pro
                end

                "Calculate property record for dynamic simulation properties using d and T as dynamic states"
                function air_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns pro::Modelica.Media.Common.ThermoFluidSpecial.ThermoProperties_dT,   # Property record for dynamic simulation
                    )
                    p::Float64,                      # Pressure
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    f = Air_Utilities.Basic.Helmholtz(d, T)
                                        pro = Modelica.Media.Common.ThermoFluidSpecial.helmholtzToProps_dT(f)
                    
                    return pro
                end

                "Calculate property record for dynamic simulation properties using p and T as dynamic states"
                function air_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns pro::Modelica.Media.Common.ThermoFluidSpecial.ThermoProperties_pT,   # Property record for dynamic simulation
                    )
                    d::Float64,                      # Density
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and dervatives w.r.t. delta and tau
                    d = Modelica.Media.Air.ReferenceAir.Air_Utilities.Inverses.dofpT(p=p, T=T, delp=1e-7)
                                        f = Air_Utilities.Basic.Helmholtz(d, T)
                                        pro = Modelica.Media.Common.ThermoFluidSpecial.helmholtzToProps_pT(f)
                    
                    return pro
                end
            end
        end
    end

    "Air: Moist air model (190 ... 647 K)"
    module MoistAir
        #=
        @extends Interfaces.PartialCondensingGases(
            mediumName="Moist air", 
            substanceNames=["water","air"], 
            reducedX=true, 
            singleState=false, 
            reference_X=[0.01,0.99], 
            fluidConstants=[IdealGases.Common.FluidData.H2O,IdealGases.Common.FluidData.N2], 
            Temperature(
            min=190, 
            max=647))
        =#
        const Water = Missing            # Index of water (in substanceNames, massFractions X, etc.)
        const Air = Missing              # Index of air (in substanceNames, massFractions X, etc.)
        const k_mair = Missing           # Ratio of molar weights
        const dryair = IdealGases.Common.DataRecord=IdealGases.Common.SingleGasesData.Air
        const steam = IdealGases.Common.DataRecord=IdealGases.Common.SingleGasesData.H2O
        const MMX = Missing              # Molar masses of components
        mutable struct ThermodynamicState
            #= @extends ThermodynamicState
             =#
        end
        @model BaseProperties begin
            #= @extends BaseProperties
            (
            T(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            p(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            Xi(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            standardOrderComponents=true) =#
            x_water                        = MassFraction(
                info="Mass of total water/mass of dry air")
            phi                            = Real(
                info="Relative humidity")
            X_liquid                       = MassFraction(
                info="Mass fraction of liquid or solid water")
            X_steam                        = MassFraction(
                info="Mass fraction of steam water")
            X_air                          = MassFraction(
                info="Mass fraction of air")
            X_sat                          = MassFraction(
                info="Steam water mass fraction of saturation boundary in kg_water/kg_moistair")
            x_sat                          = MassFraction(
                info="Steam water mass content of saturation boundary in kg_water/kg_dryair")
            p_steam_sat                    = AbsolutePressure(
                info="partial saturation pressure of steam")
            @inherits assert, T, String, mediumName, MM, Xi, Water, MMX, Air, p_steam_sat, min, saturationPressure, p, X_sat, k_mair, max, Constants, X_liquid, X_steam, X_air, h, specificEnthalpy_pTX, R, dryair, steam, u, d, state, X, x_sat, x_water, phi
            @equations begin
                        assert(T>=190 && T<=647, "
            Temperature T is not in the allowed range
            190.0 K <= (T ="+String(T)+" K) <= 647.0 K
            required from medium model \""+mediumName+"\".")
                        MM = 1/
                (Xi[Water]/MMX[Water]+(1.0-Xi[Water])/MMX[Air])
                        p_steam_sat = min(saturationPressure(T), 0.999*p)
                        X_sat = min(p_steam_sat*k_mair/max(100*Constants.eps, p-p_steam_sat)*(1-Xi[Water]), 1.0)
                        X_liquid = max(Xi[Water]-X_sat, 0.0)
                        X_steam = Xi[Water]-X_liquid
                        X_air = 1-Xi[Water]
                        h = specificEnthalpy_pTX(p, T, Xi)
                        R = dryair.R*(X_air/(1-X_liquid))+steam.R*X_steam/(1-X_liquid)
                        u = h-R*T
                        d = p/(R*T)
                        state.p = p
                        state.T = T
                        state.X = X
                        x_sat = k_mair*p_steam_sat/max(100*Constants.eps, p-p_steam_sat)
                        x_water = Xi[Water]/max(X_air, 100*Constants.eps)
                        phi = p/p_steam_sat*Xi[Water]/(Xi[Water]+k_mair*X_air)
                end
            
        end

        "Return thermodynamic state as function of pressure p, temperature T and composition X"
        function setState_pTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            T::Temperature,                  # Temperature
            X::Array{MassFraction, 1}=reference_X,   # Mass fractions
            # returns state::ThermodynamicState,   # Thermodynamic state
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=T, X=X) else ThermodynamicState(p=p, T=T, X=cat(1, X, [1-sum(X)])) end
            
            return state
        end

        "Return thermodynamic state as function of pressure p, specific enthalpy h and composition X"
        function setState_phX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            h::SpecificEnthalpy,             # Specific enthalpy
            X::Array{MassFraction, 1}=reference_X,   # Mass fractions
            # returns state::ThermodynamicState,   # Thermodynamic state
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=T_phX(p, h, X), X=X) else ThermodynamicState(p=p, T=T_phX(p, h, X), X=cat(1, X, [1-sum(X)])) end
            
            return state
        end

        "Return thermodynamic state as function of density d, temperature T and composition X"
        function setState_dTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            d::Density,                      # Density
            T::Temperature,                  # Temperature
            X::Array{MassFraction, 1}=reference_X,   # Mass fractions
            # returns state::ThermodynamicState,   # Thermodynamic state
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=d*([steam.R,dryair.R]*X)*T, T=T, X=X) else ThermodynamicState(p=d*
                ([steam.R,dryair.R]*cat(1, X, [1-sum(X)]))*T, T=T, X=cat(1, X, [1-sum(X)])) end
            
            return state
        end
        function setSmoothState(
            #= @extends setSmoothState
             =#
            )
            state = ThermodynamicState(p=Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T=Media.Common.smoothStep(x, state_a.T, state_b.T, x_small), X=Media.Common.smoothStep(x, state_a.X, state_b.X, x_small))
            
            return 
        end

        "Return absolute humidity per unit mass of moist air at saturation as a function of the thermodynamic state record"
        function Xsaturation(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns X_sat::MassFraction,   # Steam mass fraction of sat. boundary
            )
            X_sat = k_mair/
                (state.p/min(saturationPressure(state.T), 0.999*state.p)-1+k_mair)
            
            return X_sat
        end

        "Return absolute humidity per unit mass of dry air at saturation as a function of the thermodynamic state record"
        function xsaturation(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns x_sat::MassFraction,   # Absolute humidity per unit mass of dry air
            )
            x_sat = k_mair*saturationPressure(state.T)/max(100*Constants.eps, state.p-saturationPressure(state.T))
            
            return x_sat
        end

        "Return absolute humidity per unit mass of dry air at saturation as a function of pressure p and temperature T"
        function xsaturation_pT(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            T::Float64,                      # Temperature
            # returns x_sat::MassFraction,   # Absolute humidity per unit mass of dry air
            )
            x_sat = k_mair*saturationPressure(T)/max(100*Constants.eps, p-saturationPressure(T))
            
            return x_sat
        end

        "Return steam mass fraction as a function of relative humidity phi and temperature T"
        function massFraction_pTphi(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            T::Temperature,                  # Temperature
            phi::Float64,                    # Relative humidity (0 ... 1.0)
            # returns X_steam::MassFraction,   # Absolute humidity, steam mass fraction
            )
            const k::Float64 = Missing,      # Ratio of molar masses
            psat::AbsolutePressure=saturationPressure(T),   # Saturation pressure
            X_steam = phi*k/(k*phi+p/psat-phi)
            
            return X_steam
        end

        "Return relative humidity as a function of pressure p, temperature T and composition X"
        function relativeHumidity_pTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Composition
            # returns phi::Float64,          # Relative humidity
            )
            p_steam_sat::Float64,            # Saturation pressure
            X_air::Float64,                  # Dry air mass fraction
            p_steam_sat = min(saturationPressure(T), 0.999*p)
                        X_air = 1-X[Water]
                        phi = max(0.0, min(1.0, p/p_steam_sat*X[Water]/(X[Water]+k_mair*X_air)))
            
            return phi
        end

        "Return relative humidity as a function of the thermodynamic state record"
        function relativeHumidity(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state
            # returns phi::Float64,          # Relative humidity
            )
            phi = relativeHumidity_pTX(state.p, state.T, state.X)
            
            return phi
        end
        function gasConstant(
            #= @extends gasConstant
             =#
            )
            R = dryair.R*(1-state.X[Water])+steam.R*state.X[Water]
            
            return 
        end

        "Return ideal gas constant as a function from composition X"
        function gasConstant_X(
            #=
            @extends Modelica.Icons.Function()
            =#
            X::Array{Float64, 1},            # Gas phase composition
            # returns R::Float64,            # Ideal gas constant
            )
            R = dryair.R*(1-X[Water])+steam.R*X[Water]
            
            return R
        end

        "Return saturation pressure of water as a function of temperature T in the range of 273.16 to 647.096 K"
        function saturationPressureLiquid(
            #=
            @extends Modelica.Icons.Function()
            =#
            Tsat::Float64,                   # Saturation temperature
            # returns psat::Float64,         # Saturation pressure
            )
            Tcritical::Float64,              # Critical temperature
            pcritical::Float64,              # Critical pressure
            r1::Float64,                     # Common subexpression
            a::Array{Float64, 1},            # Coefficients a[:]
            n::Array{Float64, 1},            # Coefficients n[:]
            psat = exp(
                (
                (a[1]*r1^n[1]+a[2]*r1^n[2]+a[3]*r1^n[3]+a[4]*r1^n[4]+a[5]*r1^n[5]+a[6]*r1^n[6])*Tcritical)/Tsat)*pcritical
            
            return psat
        end

        "Derivative function for 'saturationPressureLiquid'"
        function saturationPressureLiquid_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            Tsat::Float64,                   # Saturation temperature
            dTsat::Float64,                  # Saturation temperature derivative
            # returns psat_der::Float64,     # Saturation pressure derivative
            )
            Tcritical::Float64,              # Critical temperature
            pcritical::Float64,              # Critical pressure
            r1::Float64,                     # Common subexpression 1
            r1_der::Float64,                 # Derivative of common subexpression 1
            a::Array{Float64, 1},            # Coefficients a[:]
            n::Array{Float64, 1},            # Coefficients n[:]
            r2::Float64,                     # Common subexpression 2
            psat_der = exp((r2*Tcritical)/Tsat)*pcritical*
                (
                (a[1]*(r1^(n[1]-1)*n[1]*r1_der)+a[2]*(r1^(n[2]-1)*n[2]*r1_der)+a[3]*(r1^(n[3]-1)*n[3]*r1_der)+a[4]*(r1^(n[4]-1)*n[4]*r1_der)+a[5]*(r1^(n[5]-1)*n[5]*r1_der)+a[6]*(r1^(n[6]-1)*n[6]*r1_der))*Tcritical/Tsat-r2*Tcritical*dTsat/Tsat^2)
            
            return psat_der
        end

        "Return sublimation pressure of water as a function of temperature T between 190 and 273.16 K"
        function sublimationPressureIce(
            #=
            @extends Modelica.Icons.Function()
            =#
            Tsat::Float64,                   # Sublimation temperature
            # returns psat::Float64,         # Sublimation pressure
            )
            Ttriple::Float64,                # Triple point temperature
            ptriple::Float64,                # Triple point pressure
            r1::Float64,                     # Common subexpression
            a::Array{Float64, 1},            # Coefficients a[:]
            n::Array{Float64, 1},            # Coefficients n[:]
            psat = exp(a[1]-a[1]*r1^n[1]+a[2]-a[2]*r1^n[2])*ptriple
            
            return psat
        end

        "Derivative function for 'sublimationPressureIce'"
        function sublimationPressureIce_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            Tsat::Float64,                   # Sublimation temperature
            dTsat::Float64,                  # Sublimation temperature derivative
            # returns psat_der::Float64,     # Sublimation pressure derivative
            )
            Ttriple::Float64,                # Triple point temperature
            ptriple::Float64,                # Triple point pressure
            r1::Float64,                     # Common subexpression 1
            r1_der::Float64,                 # Derivative of common subexpression 1
            a::Array{Float64, 1},            # Coefficients a[:]
            n::Array{Float64, 1},            # Coefficients n[:]
            psat_der = exp(a[1]-a[1]*r1^n[1]+a[2]-a[2]*r1^n[2])*ptriple*
                (-(a[1]*(r1^(n[1]-1)*n[1]*r1_der))-(a[2]*(r1^(n[2]-1)*n[2]*r1_der)))
            
            return psat_der
        end
        function saturationPressure(
            #= @extends saturationPressure
             =#
            )
            psat = Utilities.spliceFunction(saturationPressureLiquid(Tsat), sublimationPressureIce(Tsat), Tsat-273.16, 1.0)
            
            return 
        end

        "Derivative function for 'saturationPressure'"
        function saturationPressure_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            Tsat::Temperature,               # Saturation temperature
            dTsat::Float64,                  # Time derivative of saturation temperature
            # returns psat_der::Float64,     # Saturation pressure
            )
            psat_der = Utilities.spliceFunction_der(saturationPressureLiquid(Tsat), sublimationPressureIce(Tsat), Tsat-273.16, 1.0, saturationPressureLiquid_der(Tsat=Tsat, dTsat=dTsat), sublimationPressureIce_der(Tsat=Tsat, dTsat=dTsat), dTsat, 0)
            
            return psat_der
        end

        "Return saturation temperature of water as a function of (partial) pressure p"
        function saturationTemperature(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T_min::Float64,                  # Lower boundary of solution
            T_max::Float64,                  # Upper boundary of solution
            # returns T::Float64,            # Saturation temperature
            )
#=

            
            module Internal
                #=
                @extends Modelica.Media.Common.OneNonLinearEquation()
                =#
                mutable struct f_nonlinear_Data
                    #= @extends f_nonlinear_Data
                     =#
                end
                function f_nonlinear(
                    #= @extends f_nonlinear
                     =#
                    )
                    y = saturationPressure(x)
                    
                    return 
                end
                function solve(
                    #= @extends solve
                     =#
                    )
                    return 
                end
            end
=#
            T = Internal.solve(p, T_min, T_max, f_nonlinear_data=Internal.f_nonlinear_Data())
            
            return T
        end
        function enthalpyOfVaporization(
            #= @extends enthalpyOfVaporization
             =#
            )
            Tcritical::Float64,              # Critical temperature
            dcritical::Float64,              # Critical density
            pcritical::Float64,              # Critical pressure
            n::Array{Float64, 1},            # Powers in equation (1)
            a::Array{Float64, 1},            # Coefficients in equation (1) of [1]
            m::Array{Float64, 1},            # Powers in equation (2)
            b::Array{Float64, 1},            # Coefficients in equation (2) of [1]
            o::Array{Float64, 1},            # Powers in equation (3)
            c::Array{Float64, 1},            # Coefficients in equation (3) of [1]
            tau::Float64,                    # Temperature expression
            r1::Float64,                     # Expression 1
            r2::Float64,                     # Expression 2
            dp::Float64,                     # Density of saturated liquid
            dpp::Float64,                    # Density of saturated vapor
            r0 = -
                (((dp-dpp)*exp(r1)*pcritical*(r2+r1*tau))/(dp*dpp*tau))
            
            return 
        end

        "Return specific heat capacity of water (liquid only) as a function of temperature T"
        function HeatCapacityOfWater(
            #=
            @extends Modelica.Icons.Function()
            =#
            T::Temperature,                  # Temperature
            # returns cp_fl::SpecificHeatCapacity,   # Specific heat capacity of liquid
            )
            cp_fl = 1e3*
                (4.2166-(T-273.15)*
                (0.0033166+(T-273.15)*
                (0.00010295-(T-273.15)*(1.3819e-6+(T-273.15)*7.3221e-9))))
            
            return cp_fl
        end
        function enthalpyOfLiquid(
            #= @extends enthalpyOfLiquid
             =#
            )
            h = (T-273.15)*1e3*
                (4.2166-0.5*(T-273.15)*
                (0.0033166+0.333333*(T-273.15)*
                (0.00010295-0.25*(T-273.15)*(1.3819e-6+0.2*(T-273.15)*7.3221e-9))))
            
            return 
        end
        function enthalpyOfGas(
            #= @extends enthalpyOfGas
             =#
            )
            h = Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)*X[Water]+Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)*(1.0-X[Water])
            
            return 
        end
        function enthalpyOfCondensingGas(
            #= @extends enthalpyOfCondensingGas
             =#
            )
            h = Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)
            
            return 
        end
        function enthalpyOfNonCondensingGas(
            #= @extends enthalpyOfNonCondensingGas
             =#
            )
            h = Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)
            
            return 
        end

        "Computes specific enthalpy of water (solid/liquid) near atmospheric pressure from temperature T"
        function enthalpyOfWater(
            #=
            @extends Modelica.Icons.Function()
            =#
            T::Float64,                      # Temperature
            # returns h::Float64,            # Specific enthalpy of water
            )
            h = Utilities.spliceFunction(4200*(T-273.15), 2050*(T-273.15)-333000, T-273.16, 0.1)
            
            return h
        end

        "Derivative function of enthalpyOfWater"
        function enthalpyOfWater_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            T::Float64,                      # Temperature
            dT::Float64,                     # Time derivative of temperature
            # returns dh::Float64,           # Time derivative of specific enthalpy
            )
            dh = Utilities.spliceFunction_der(4200*(T-273.15), 2050*(T-273.15)-333000, T-273.16, 0.1, 4200*dT, 2050*dT, dT, 0)
            
            return dh
        end
        function pressure(
            #= @extends pressure
             =#
            )
            p = state.p
            
            return 
        end
        function temperature(
            #= @extends temperature
             =#
            )
            T = state.T
            
            return 
        end

        "Return temperature as a function of pressure p, specific enthalpy h and composition X"
        function T_phX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            h::SpecificEnthalpy,             # Specific enthalpy
            X::Array{MassFraction, 1},       # Mass fractions of composition
            # returns T::Temperature,        # Temperature
            )
#=

            "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
            module Internal
                #=
                @extends Modelica.Media.Common.OneNonLinearEquation()
                =#
                mutable struct f_nonlinear_Data
                    #= @extends f_nonlinear_Data
                     =#
                    #=
                    @extends Modelica.Media.IdealGases.Common.DataRecord()
                    =#
                end
                function f_nonlinear(
                    #= @extends f_nonlinear
                     =#
                    )
                    y = h_pTX(p, x, X)
                    
                    return 
                end
                function solve(
                    #= @extends solve
                     =#
                    )
                    return 
                end
            end
=#
            T = Internal.solve(h, 190, 647, p, X[1:nXi], steam)
            
            return T
        end
        function density(
            #= @extends density
             =#
            )
            d = state.p/(gasConstant(state)*state.T)
            
            return 
        end
        function specificEnthalpy(
            #= @extends specificEnthalpy
             =#
            )
            h = h_pTX(state.p, state.T, state.X)
            
            return 
        end

        "Return specific enthalpy of moist air as a function of pressure p, temperature T and composition X"
        function h_pTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            # returns h::Float64,            # Specific enthalpy at p, T, X
            )
            p_steam_sat::Float64,            # partial saturation pressure of steam
            X_sat::Float64,                  # Absolute humidity per unit mass of moist air
            X_liquid::Float64,               # Mass fraction of liquid water
            X_steam::Float64,                # Mass fraction of steam water
            X_air::Float64,                  # Mass fraction of air
            p_steam_sat = saturationPressure(T)
                        X_sat = min(p_steam_sat*k_mair/max(100*Constants.eps, p-p_steam_sat)*(1-X[Water]), 1.0)
                        X_liquid = max(X[Water]-X_sat, 0.0)
                        X_steam = X[Water]-X_liquid
                        X_air = 1-X[Water]
                        h = [Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5),Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)]*[X_steam,X_air]+enthalpyOfWater(T)*X_liquid
            
            return h
        end

        "Derivative function of h_pTX"
        function h_pTX_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            dp::Float64,                     # Pressure derivative
            dT::Float64,                     # Temperature derivative
            dX::Array{Float64, 1},           # Composition derivative
            # returns h_der::Float64,        # Time derivative of specific enthalpy
            )
            p_steam_sat::Float64,            # partial saturation pressure of steam
            X_sat::Float64,                  # Absolute humidity per unit mass of moist air
            X_liquid::Float64,               # Mass fraction of liquid water
            X_steam::Float64,                # Mass fraction of steam water
            X_air::Float64,                  # Mass fraction of air
            x_sat::Float64,                  # Absolute humidity per unit mass of dry air at saturation
            dX_steam::Float64,               # Time derivative of steam mass fraction
            dX_air::Float64,                 # Time derivative of dry air mass fraction
            dX_liq::Float64,                 # Time derivative of liquid/solid water mass fraction
            dps::Float64,                    # Time derivative of saturation pressure
            dx_sat::Float64,                 # Time derivative of absolute humidity per unit mass of dry air
            p_steam_sat = saturationPressure(T)
                        x_sat = p_steam_sat*k_mair/max(100*Modelica.Constants.eps, p-p_steam_sat)
                        X_sat = min(x_sat*(1-X[Water]), 1.0)
                        X_liquid = Utilities.smoothMax(X[Water]-X_sat, 0.0, 1e-5)
                        X_steam = X[Water]-X_liquid
                        X_air = 1-X[Water]
                        dX_air = -dX[Water]
                        dps = saturationPressure_der(Tsat=T, dTsat=dT)
                        dx_sat = k_mair*
                (dps*(p-p_steam_sat)-p_steam_sat*(dp-dps))/(p-p_steam_sat)/(p-p_steam_sat)
                        dX_liq = Utilities.smoothMax_der(X[Water]-X_sat, 0.0, 1e-5, (1+x_sat)*dX[Water]-(1-X[Water])*dx_sat, 0, 0)
                        dX_steam = dX[Water]-dX_liq
                        h_der = X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5, dT=dT)+dX_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)+X_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684, dT=dT)+dX_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)+X_liquid*enthalpyOfWater_der(T=T, dT=dT)+dX_liq*enthalpyOfWater(T)
            
            return h_der
        end
        function isentropicExponent(
            #= @extends isentropicExponent
             =#
            )
            gamma = specificHeatCapacityCp(state)/specificHeatCapacityCv(state)
            
            return 
        end

        "Approximate calculation of h_is from upstream properties, downstream pressure, gas part only"
        function isentropicEnthalpyApproximation(
            #=
            @extends Modelica.Icons.Function()
            =#
            p2::AbsolutePressure,            # Downstream pressure
            state::ThermodynamicState,       # Thermodynamic state at upstream location
            # returns h_is::SpecificEnthalpy,   # Isentropic enthalpy
            )
            h::SpecificEnthalpy,             # Specific enthalpy at upstream location
            gamma::IsentropicExponent=isentropicExponent(state),   # Isentropic exponent
            X::Array{MassFraction, 1},       # Complete X-vector
            X = state.X
                        h = [Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=state.T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5),Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=state.T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)]*X
                        h_is = h+gamma/(gamma-1.0)*(state.T*gasConstant(state))*((p2/state.p)^((gamma-1)/gamma)-1.0)
            
            return h_is
        end
        function specificInternalEnergy(
            #= @extends specificInternalEnergy
             =#
            #=
            @extends Modelica.Icons.Function()
            =#
            # returns u::Float64,            # Specific internal energy
            )
            u = specificInternalEnergy_pTX(state.p, state.T, state.X)
            
            return u
        end

        "Return specific internal energy of moist air as a function of pressure p, temperature T and composition X"
        function specificInternalEnergy_pTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            # returns u::Float64,            # Specific internal energy
            )
            p_steam_sat::Float64,            # partial saturation pressure of steam
            X_liquid::Float64,               # Mass fraction of liquid water
            X_steam::Float64,                # Mass fraction of steam water
            X_air::Float64,                  # Mass fraction of air
            X_sat::Float64,                  # Absolute humidity per unit mass of moist air
            R_gas::Float64,                  # Ideal gas constant
            p_steam_sat = saturationPressure(T)
                        X_sat = min(p_steam_sat*k_mair/max(100*Constants.eps, p-p_steam_sat)*(1-X[Water]), 1.0)
                        X_liquid = max(X[Water]-X_sat, 0.0)
                        X_steam = X[Water]-X_liquid
                        X_air = 1-X[Water]
                        R_gas = dryair.R*X_air/(1-X_liquid)+steam.R*X_steam/(1-X_liquid)
                        u = X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)+X_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)+enthalpyOfWater(T)*X_liquid-R_gas*T
            
            return u
        end

        "Derivative function for specificInternalEnergy_pTX"
        function specificInternalEnergy_pTX_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            dp::Float64,                     # Pressure derivative
            dT::Float64,                     # Temperature derivative
            dX::Array{Float64, 1},           # Mass fraction derivatives
            # returns u_der::Float64,        # Specific internal energy derivative
            )
            p_steam_sat::Float64,            # partial saturation pressure of steam
            X_liquid::Float64,               # Mass fraction of liquid water
            X_steam::Float64,                # Mass fraction of steam water
            X_air::Float64,                  # Mass fraction of air
            X_sat::Float64,                  # Absolute humidity per unit mass of moist air
            R_gas::Float64,                  # Ideal gas constant
            x_sat::Float64,                  # Absolute humidity per unit mass of dry air at saturation
            dX_steam::Float64,               # Time derivative of steam mass fraction
            dX_air::Float64,                 # Time derivative of dry air mass fraction
            dX_liq::Float64,                 # Time derivative of liquid/solid water mass fraction
            dps::Float64,                    # Time derivative of saturation pressure
            dx_sat::Float64,                 # Time derivative of absolute humidity per unit mass of dry air
            dR_gas::Float64,                 # Time derivative of ideal gas constant
            p_steam_sat = saturationPressure(T)
                        x_sat = p_steam_sat*k_mair/max(100*Modelica.Constants.eps, p-p_steam_sat)
                        X_sat = min(x_sat*(1-X[Water]), 1.0)
                        X_liquid = Utilities.spliceFunction(X[Water]-X_sat, 0.0, X[Water]-X_sat, 1e-6)
                        X_steam = X[Water]-X_liquid
                        X_air = 1-X[Water]
                        R_gas = steam.R*X_steam/(1-X_liquid)+dryair.R*X_air/(1-X_liquid)
                        dX_air = -dX[Water]
                        dps = saturationPressure_der(Tsat=T, dTsat=dT)
                        dx_sat = k_mair*
                (dps*(p-p_steam_sat)-p_steam_sat*(dp-dps))/(p-p_steam_sat)/(p-p_steam_sat)
                        dX_liq = Utilities.spliceFunction_der(X[Water]-X_sat, 0.0, X[Water]-X_sat, 1e-6, (1+x_sat)*dX[Water]-(1-X[Water])*dx_sat, 0.0, (1+x_sat)*dX[Water]-(1-X[Water])*dx_sat, 0.0)
                        dX_steam = dX[Water]-dX_liq
                        dR_gas = 
                (steam.R*(dX_steam*(1-X_liquid)+dX_liq*X_steam)+dryair.R*(dX_air*(1-X_liquid)+dX_liq*X_air))/(1-X_liquid)/(1-X_liquid)
                        u_der = X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5, dT=dT)+dX_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=steam, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5)+X_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684, dT=dT)+dX_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)+X_liquid*enthalpyOfWater_der(T=T, dT=dT)+dX_liq*enthalpyOfWater(T)-dR_gas*T-R_gas*dT
            
            return u_der
        end
        function specificEntropy(
            #= @extends specificEntropy
             =#
            )
            s = s_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificGibbsEnergy(
            #= @extends specificGibbsEnergy
             =#
            #=
            @extends Modelica.Icons.Function()
            =#
            )
            g = h_pTX(state.p, state.T, state.X)-state.T*specificEntropy(state)
            
            return 
        end
        function specificHelmholtzEnergy(
            #= @extends specificHelmholtzEnergy
             =#
            #=
            @extends Modelica.Icons.Function()
            =#
            )
            f = h_pTX(state.p, state.T, state.X)-gasConstant(state)*state.T-state.T*specificEntropy(state)
            
            return 
        end
        function specificHeatCapacityCp(
            #= @extends specificHeatCapacityCp
             =#
            )
            dT::Float64,                  
            cp = h_pTX_der(state.p, state.T, state.X, 0.0, 1.0, zeros(size(state.X, 1)))*dT
            
            return 
        end
        function specificHeatCapacityCv(
            #= @extends specificHeatCapacityCv
             =#
            )
            cv = Modelica.Media.IdealGases.Common.Functions.cp_Tlow(dryair, state.T)*(1-state.X[Water])+Modelica.Media.IdealGases.Common.Functions.cp_Tlow(steam, state.T)*state.X[Water]-gasConstant(state)
            
            return 
        end
        function dynamicViscosity(
            #= @extends dynamicViscosity
             =#
            )
            eta = 1e-6*Polynomials_Temp.evaluateWithRange([9.7391102886305869E-15,-3.1353724870333906E-11,4.3004876595642225E-08,-3.8228016291758240E-05,5.0427874367180762E-02,1.7239260139242528E+01], Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T))
            
            return 
        end
        function thermalConductivity(
            #= @extends thermalConductivity
             =#
            )
            lambda = 1e-3*Polynomials_Temp.evaluateWithRange([6.5691470817717812E-15,-3.4025961923050509E-11,5.3279284846303157E-08,-4.5340839289219472E-05,7.6129675309037664E-02,2.4169481088097051E+01], Cv.to_degC(123.15), Cv.to_degC(1273.15), Cv.to_degC(state.T))
            
            return 
        end
        function velocityOfSound(
            #= @extends velocityOfSound
             =#
            )
            a = sqrt(isentropicExponent(state)*gasConstant(state)*temperature(state))
            
            return 
        end
        function isobaricExpansionCoefficient(
            #= @extends isobaricExpansionCoefficient
             =#
            )
            beta = 1/temperature(state)
            
            return 
        end
        function isothermalCompressibility(
            #= @extends isothermalCompressibility
             =#
            )
            kappa = 1/pressure(state)
            
            return 
        end
        function density_derp_h(
            #= @extends density_derp_h
             =#
            )
            ddph = 1/(gasConstant(state)*temperature(state))
            
            return 
        end
        function density_derh_p(
            #= @extends density_derh_p
             =#
            )
            ddhp = -density(state)/
                (specificHeatCapacityCp(state)*temperature(state))
            
            return 
        end
        function density_derp_T(
            #= @extends density_derp_T
             =#
            )
            ddpT = 1/(gasConstant(state)*temperature(state))
            
            return 
        end
        function density_derT_p(
            #= @extends density_derT_p
             =#
            )
            ddTp = -density(state)/temperature(state)
            
            return 
        end
        function density_derX(
            #= @extends density_derX
             =#
            )
            dddX[Water] = pressure(state)*(steam.R-dryair.R)/
                ((steam.R-dryair.R)*state.X[Water]*temperature(state)+dryair.R*temperature(state))^2
                        dddX[Air] = pressure(state)*(dryair.R-steam.R)/
                ((dryair.R-steam.R)*state.X[Air]*temperature(state)+steam.R*temperature(state))^2
            
            return 
        end
        function molarMass(
            #= @extends molarMass
             =#
            )
            MM = Modelica.Media.Air.MoistAir.gasConstant(state)/Modelica.Constants.R
            
            return 
        end

        "Return temperature as a function of pressure p, specific entropy s and composition X"
        function T_psX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Pressure
            s::SpecificEntropy,              # Specific entropy
            X::Array{MassFraction, 1},       # Mass fractions of composition
            # returns T::Temperature,        # Temperature
            )
#=

            "Solve s(data,T) for T with given s"
            module Internal
                #=
                @extends Modelica.Media.Common.OneNonLinearEquation()
                =#
                mutable struct f_nonlinear_Data
                    #= @extends f_nonlinear_Data
                     =#
                    #=
                    @extends Modelica.Media.IdealGases.Common.DataRecord()
                    =#
                end
                function f_nonlinear(
                    #= @extends f_nonlinear
                     =#
                    )
                    y = s_pTX(p, x, X)
                    
                    return 
                end
                function solve(
                    #= @extends solve
                     =#
                    )
                    return 
                end
            end
=#
            T = Internal.solve(s, 190, 647, p, X[1:nX], steam)
            
            return T
        end
        function setState_psX(
            #= @extends setState_psX
             =#
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=T_psX(p, s, X), X=X) else ThermodynamicState(p=p, T=T_psX(p, s, X), X=cat(1, X, [1-sum(X)])) end
            
            return 
        end

        "Return specific entropy of moist air as a function of pressure p, temperature T and composition X (only valid for phi<1)"
        function s_pTX(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            # returns s::Float64,            # Specific entropy at p, T, X
            )
            Y::Array{MoleFraction, 1}=massToMoleFractions(X, [steam.MM,dryair.MM]),   # Molar fraction
            s = Modelica.Media.IdealGases.Common.Functions.s0_Tlow(dryair, T)*(1-X[Water])+Modelica.Media.IdealGases.Common.Functions.s0_Tlow(steam, T)*X[Water]-Modelica.Constants.R*
                (Utilities.smoothMax(X[Water]/MMX[Water]*Modelica.Math.log(max(Y[Water], Modelica.Constants.eps)*p/reference_p), 0.0, 1e-9)-Utilities.smoothMax((1-X[Water])/MMX[Air]*Modelica.Math.log(max(Y[Air], Modelica.Constants.eps)*p/reference_p), 0.0, 1e-9))
            
            return s
        end

        "Return specific entropy of moist air as a function of pressure p, temperature T and composition X (only valid for phi<1)"
        function s_pTX_der(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::Float64,                      # Pressure
            T::Float64,                      # Temperature
            X::Array{Float64, 1},            # Mass fractions of moist air
            dp::Float64,                     # Derivative of pressure
            dT::Float64,                     # Derivative of temperature
            dX::Array{Float64, 1},           # Derivative of mass fractions
            # returns ds::Float64,           # Specific entropy at p, T, X
            )
            Y::Array{MoleFraction, 1}=massToMoleFractions(X, [steam.MM,dryair.MM]),   # Molar fraction
            ds = Modelica.Media.IdealGases.Common.Functions.s0_Tlow_der(dryair, T, dT)*(1-X[Water])+Modelica.Media.IdealGases.Common.Functions.s0_Tlow_der(steam, T, dT)*X[Water]+Modelica.Media.IdealGases.Common.Functions.s0_Tlow(dryair, T)*dX[Air]+Modelica.Media.IdealGases.Common.Functions.s0_Tlow(steam, T)*dX[Water]-Modelica.Constants.R*
                (1/MMX[Water]*Utilities.smoothMax_der(X[Water]*Modelica.Math.log(max(Y[Water], Modelica.Constants.eps)*p/reference_p), 0.0, 1e-9, 
                (Modelica.Math.log(max(Y[Water], Modelica.Constants.eps)*p/reference_p)+
                (X[Water]/Y[Water]*
                (X[Air]*MMX[Water]/(X[Air]*MMX[Water]+X[Water]*MMX[Air])^2)))*dX[Water]+X[Water]*reference_p/p*dp, 0, 0)-1/MMX[Air]*Utilities.smoothMax_der((1-X[Water])*Modelica.Math.log(max(Y[Air], Modelica.Constants.eps)*p/reference_p), 0.0, 1e-9, 
                (Modelica.Math.log(max(Y[Air], Modelica.Constants.eps)*p/reference_p)+
                (X[Air]/Y[Air]*
                (X[Water]*MMX[Air]/(X[Air]*MMX[Water]+X[Water]*MMX[Air])^2)))*dX[Air]+X[Air]*reference_p/p*dp, 0, 0))
            
            return ds
        end
        function isentropicEnthalpy(
            #= @extends isentropicEnthalpy
             =#
            #=
            @extends Modelica.Icons.Function()
            =#
            )
            h_is = Modelica.Media.Air.MoistAir.h_pTX(p_downstream, Modelica.Media.Air.MoistAir.T_psX(p_downstream, Modelica.Media.Air.MoistAir.specificEntropy(refState), refState.X), refState.X)
            
            return 
        end

        "Utility functions"
        module Utilities
            #=
            @extends Modelica.Icons.UtilitiesPackage()
            =#

            "Spline interpolation of two functions"
            function spliceFunction(
                #=
                @extends Modelica.Icons.Function()
                =#
                pos::Float64,                    # Returned value for x-deltax >= 0
                neg::Float64,                    # Returned value for x+deltax <= 0
                x::Float64,                      # Function argument
                deltax::Float64,                 # Region around x with spline interpolation
                # returns out::Float64,       
                )
                scaledX::Float64,             
                scaledX1::Float64,            
                y::Float64,                   
                scaledX1 = x/deltax
                                scaledX = scaledX1*Modelica.Math.asin(1)
                                if scaledX1<=-0.999999999
                    y = 0
                elseif scaledX1>=0.999999999
                    y = 1
                else
                    y = 
                    (Modelica.Math.tanh(Modelica.Math.tan(scaledX))+1)/2
                end
                                out = pos*y+(1-y)*neg
                
                return out
            end

            "Derivative of spliceFunction"
            function spliceFunction_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                pos::Float64,                 
                neg::Float64,                 
                x::Float64,                   
                deltax::Float64,              
                dpos::Float64,                
                dneg::Float64,                
                dx::Float64,                  
                ddeltax::Float64,             
                # returns out::Float64,       
                )
                scaledX::Float64,             
                scaledX1::Float64,            
                dscaledX1::Float64,           
                y::Float64,                   
                scaledX1 = x/deltax
                                scaledX = scaledX1*Modelica.Math.asin(1)
                                dscaledX1 = (dx-scaledX1*ddeltax)/deltax
                                if scaledX1<=-0.99999999999
                    y = 0
                elseif scaledX1>=0.9999999999
                    y = 1
                else
                    y = 
                    (Modelica.Math.tanh(Modelica.Math.tan(scaledX))+1)/2
                end
                                out = dpos*y+(1-y)*dneg
                                if (abs(scaledX1)<1)
                    out = out+(pos-neg)*dscaledX1*Modelica.Math.asin(1)/2/
                    (Modelica.Math.cosh(Modelica.Math.tan(scaledX))*Modelica.Math.cos(scaledX))^2
                end
                
                return out
            end

            
            function smoothMax(
                #=
                @extends Modelica.Icons.Function()
                =#
                x1::Float64,                     # First argument of smooth max operator
                x2::Float64,                     # Second argument of smooth max operator
                dx::Float64,                     # Approximate difference between x1 and x2, below which regularization starts
                # returns y::Float64,            # Result of smooth max operator
                )
                y = max(x1, x2)+Math.log((exp((4/dx)*(x1-max(x1, x2))))+(exp((4/dx)*(x2-max(x1, x2)))))/(4/dx)
                
                return y
            end

            
            function smoothMax_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                x1::Float64,                     # First argument of smooth max operator
                x2::Float64,                     # Second argument of smooth max operator
                dx::Float64,                     # Approximate difference between x1 and x2, below which regularization starts
                dx1::Float64,                 
                dx2::Float64,                 
                ddx::Float64,                 
                # returns dy::Float64,           # Derivative of smooth max operator
                )
                dy = (if x1>x2; dx1 else dx2 end)+0.25*
                    (
                    (
                    (4*(dx1-(if x1>x2; dx1 else dx2 end))/dx-4*(x1-max(x1, x2))*ddx/dx^2)*exp(4*(x1-max(x1, x2))/dx)+
                    (4*(dx2-(if x1>x2; dx1 else dx2 end))/dx-4*(x2-max(x1, x2))*ddx/dx^2)*exp(4*(x2-max(x1, x2))/dx))*dx/
                    (exp(4*(x1-max(x1, x2))/dx)+exp(4*(x2-max(x1, x2))/dx))+log(exp(4*(x1-max(x1, x2))/dx)+exp(4*(x2-max(x1, x2))/dx))*ddx)
                
                return dy
            end
        end
    end

    "ReferenceMoistAir: Detailed moist air model (143.15 ... 2000 K)"
    module ReferenceMoistAir
        #=
        @extends Modelica.Media.Interfaces.PartialRealCondensingGases(
            mediumName="Moist air", 
            substanceNames=["Water","Air"], 
            fixedX=false, 
            reducedX=true, 
            singleState=false, 
            reference_X=[0.01,0.99], 
            fluidConstants=[Utilities.Water95_Utilities.waterConstants,Modelica.Media.Air.ReferenceAir.airConstants], 
            ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX)
        =#
        const Water = Missing            # Index of water (in substanceNames, massFractions X, etc.)
        const Air = Missing              # Index of air (in substanceNames, massFractions X, etc.)
        const useEnhancementFactor = Missing   # Use the enhancement factor in the calculations
        const useDissociation = Missing   # Take dissociation into account for high temperatures
        const k_mair = Missing           # Ratio of molar weights
        const dryair = Common.FundamentalConstants=ReferenceAir.Air_Utilities.Basic.Constants
        const steam = Common.FundamentalConstants=Utilities.Water95_Utilities.Constants
        const MMX = Missing              # Molar masses of components
        mutable struct ThermodynamicState
            #= @extends ThermodynamicState
             =#
        end
        @model BaseProperties begin
            #= @extends BaseProperties
            (
            T(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            p(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            Xi(
            stateSelect=if preferredMediumStates; StateSelect.prefer else StateSelect.default end), 
            standardOrderComponents=true) =#
            x_water                        = MassFraction(
                info="Mass of total water/mass of dry air")
            phi                            = Real(
                info="Relative humidity")
            X_liquid                       = MassFraction(
                info="Mass fraction of liquid or solid water")
            X_steam                        = MassFraction(
                info="Mass fraction of steam water")
            X_air                          = MassFraction(
                info="Mass fraction of air")
            X_sat                          = MassFraction(
                info="Steam water mass fraction of saturation boundary in kg_water/kg_moistair")
            x_sat                          = MassFraction(
                info="Steam water mass content of saturation boundary in kg_water/kg_dryair")
            p_steam_sat                    = AbsolutePressure(
                info="partial saturation pressure of steam")
            @inherits assert, T, String, mediumName, MM, Xi, Water, MMX, Air, p_steam_sat, Modelica, p, X_sat, k_mair, X_liquid, X_steam, X_air, h, specificEnthalpy_pTX, R, dryair, steam, u, d, state, X, x_sat, max, Constants, x_water, phi
            @equations begin
                        assert(T>=143.15 && T<=2000, "Temperature T is not in the allowed range 143.15 K <= (T ="+String(T)+" K) <= 2000 K required from medium model \""+mediumName+"\".")
                        MM = 1/
                (Xi[Water]/MMX[Water]+(1.0-Xi[Water])/MMX[Air])
                        p_steam_sat = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p, T)
                        X_sat = k_mair/(p/p_steam_sat-1+k_mair)
                        X_liquid = Xi[Water]-X_sat
                        X_steam = Xi[Water]-X_liquid
                        X_air = 1-Xi[Water]
                        h = specificEnthalpy_pTX(p, T, Xi)
                        R = dryair.R*(X_air/(1-X_liquid))+steam.R*X_steam/(1-X_liquid)
                        u = Modelica.Media.Air.ReferenceMoistAir.Utilities.u_pTX(p, T, Xi)
                        d = Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX(p, T, Xi)
                        state.p = p
                        state.T = T
                        state.X = X
                        x_sat = k_mair*p_steam_sat/max(100*Constants.eps, p-p_steam_sat)
                        x_water = Xi[Water]/max(X_air, 100*Constants.eps)
                        phi = Modelica.Media.Air.ReferenceMoistAir.Utilities.phi_pTX(p, T, Xi)
                end
            
        end
        function setState_pTX(
            #= @extends setState_pTX
             =#
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=T, X=X) else ThermodynamicState(p=p, T=T, X=cat(1, X, [1-sum(X)])) end
            
            return 
        end
        function setState_phX(
            #= @extends setState_phX
             =#
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.T_phX(p, h, X), X=X) else ThermodynamicState(p=p, T=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.T_phX(p, h, X), X=cat(1, X, [1-sum(X)])) end
            
            return 
        end
        function setState_psX(
            #= @extends setState_psX
             =#
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=p, T=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.T_psX(p, s, X), X=X) else ThermodynamicState(p=p, T=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.T_psX(p, s, X), X=cat(1, X, [1-sum(X)])) end
            
            return 
        end
        function setState_dTX(
            #= @extends setState_dTX
             =#
            )
            state = if size(X, 1)==nX; ThermodynamicState(p=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.p_dTX(d, T, X), T=T, X=X) else ThermodynamicState(p=Modelica.Media.Air.ReferenceMoistAir.Utilities.Inverses.p_dTX(d, T, X), T=T, X=cat(1, X, [1-sum(X)])) end
            
            return 
        end
        function setSmoothState(
            #= @extends setSmoothState
             =#
            )
            state = ThermodynamicState(p=Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small), T=Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small), X=Modelica.Media.Common.smoothStep(x, state_a.X, state_b.X, x_small))
            
            return 
        end

        "Return absolute humitity per unit mass of moist air at saturation as a function of the thermodynamic state record"
        function Xsaturation(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns X_sat::MassFraction,   # Steam mass fraction of sat. boundary
            )
            X::Array{MassFraction, 1},    
            X = massFractionSaturation(state)
                        X_sat = X[1]
            
            return X_sat
        end

        "Return absolute humitity per unit mass of dry air at saturation as a function of the thermodynamic state record"
        function xsaturation(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns x_sat::MassFraction,   # Absolute humidity per unit mass of dry air
            )
            x_sat = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(state.p, state.T)
                        assert(x_sat>-1, "Calculation of absolute humidity is meaningless\nfor input pressure p = "+String(state.p)+" Pa and temperature T = "+String(state.T)+" K.")
            
            return x_sat
        end
        function massFraction_pTphi(
            #= @extends massFraction_pTphi
             =#
            )
            pds::Float64,                 
            assert(phi<1.0 && phi>0, "Illegal input phi = "+String(phi)+". Relative humidity is only defined in the range\n 0 <= phi <= 1.0.")
                        pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p, T)
                        assert(pds>-1, "Calculation of mass fraction of steam is meaningless\nfor input pressure p = "+String(p)+" Pa and temperature T = "+String(T)+" K.")
                        X = [phi*k_mair/(p/pds-phi),1-phi*k_mair/(p/pds-phi)]
            
            return 
        end

        "Return mass fraction of water vapor"
        function massFractionWaterVapor(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns X::MassFraction,       # Mass fraction of water vapor
            )
            xw::Float64,                  
            xws::Float64,                 
            xw = state.X[1]/(1-state.X[1])
                        xws = Utilities.xws_pT(state.p, state.T)
                        X = if (xw<=xws); xw/(1+xw) else xws/(1+xw) end
            
            return X
        end

        "Return mass fraction of liquid and solid water"
        function massFractionWaterNonVapor(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns X::MassFraction,       # Mass fraction of water varpor
            )
            xw::Float64,                  
            xws::Float64,                 
            xw = state.X[1]/(1-state.X[1])
                        xws = Utilities.xws_pT(state.p, state.T)
                        X = if (xw<=xws); 0 else (xw-xws)/(1+xw) end
            
            return X
        end
        function massFractionSaturation(
            #= @extends massFractionSaturation
             =#
            )
            pds::AbsolutePressure,        
            pds = Utilities.pds_pT(state.p, state.T)
                        Xsat = [k_mair/(state.p/pds-1+k_mair),(state.p/pds-1)/(state.p/pds-1+k_mair)]
                        assert(Xsat[1]>-1, "Calculation of saturation mass fraction is meaningless\nfor input pressure p = "+String(state.p)+" Pa and temperature T = "+String(state.T)+" K.")
            
            return 
        end

        "Return mass fvraction at saturation boundary given pressure and saturation pressure"
        function massFractionSaturation_ppsat(
            #=
            @extends Modelica.Icons.Function()
            =#
            p::AbsolutePressure,             # Ambient pressure
            psat::AbsolutePressure,          # Saturation pressure
            # returns X::Array{MassFraction, 1},   # Mass fraction
            )
            X = [k_mair/(p/psat-1+k_mair),(p/psat-1)/(p/psat-1+k_mair)]
            
            return X
        end

        "Return mass fractions as a function of pressure, temperature and absolute humidity in kg(water)/kg(dry air)"
        function massFraction_waterContent(
            #=
            @extends Modelica.Icons.Function()
            =#
            xw::Float64,                     # Water content in kg(water)/kg(dry air)
            # returns X::Array{MassFraction, 1},   # Mass fractions
            )
            X = [xw/(1+xw),1/(1+xw)]
            
            return X
        end

        "Return water content in kg(water)/kg(dry air) given mass fractions"
        function waterContent_X(
            #=
            @extends Modelica.Icons.Function()
            =#
            X::Array{MassFraction, 1},       # Mass fractions
            # returns xw::Float64,           # Water content in kg(water)/kg(dry air)
            )
            xw = X[1]/(1-X[1])
            
            return xw
        end
        function relativeHumidity(
            #= @extends relativeHumidity
             =#
            )
            phi = Utilities.phi_pTX(state.p, state.T, state.X)
                        assert(phi>-1, "Calculation of relative humidity is meaningless\nfor input pressure p = "+String(state.p)+" Pa and temperature T = "+String(state.T)+" K.")
            
            return 
        end
        function gasConstant(
            #= @extends gasConstant
             =#
            )
            R = dryair.R*(1-state.X[Water])+steam.R*state.X[Water]
            
            return 
        end

        "Return saturation pressure of water as a function of temperature T"
        function saturationPressureLiquid(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns psat::AbsolutePressure,   # Saturation pressure
            )
            psat = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(state.T)
            
            return psat
        end

        "Return sublimation pressure of water as a function of temperature T between 223.16 and 273.16 K"
        function sublimationPressureIce(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns psat::AbsolutePressure,   # Sublimation pressure
            )
            psat = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub(state.T)
            
            return psat
        end
        function saturationPressure(
            #= @extends saturationPressure
             =#
            )
            psat = Utilities.pds_pT(state.p, state.T)
                        assert(psat>-1, "Calculation of saturation pressure is meaningless\nfor input temperature T = "+String(state.T)+" K.")
            
            return 
        end
        function saturationTemperature(
            #= @extends saturationTemperature
             =#
            )
#=

            
            function Tsat_res(
                #=
                @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                =#
                p::Float64,                      # Pressure
                )
                y = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p=p, T=u)-p
                
                return 
            end
=#
            Tsat = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
            
            return 
        end
        function enthalpyOfVaporization(
            #= @extends enthalpyOfVaporization
             =#
            )
            p_liq::AbsolutePressure,      
            p_liq = saturationPressureLiquid(state)
                        r0 = Modelica.Media.Water.IF97_Utilities.hv_p(p_liq)-Modelica.Media.Water.IF97_Utilities.hl_p(p_liq)
            
            return 
        end
        function enthalpyOfLiquid(
            #= @extends enthalpyOfLiquid
             =#
            )
            xw::Float64,                  
            xws::Float64,                 
            xw = state.X[1]/(1-state.X[1])
                        xws = Utilities.xws_pT(state.p, state.T)
                        if ((xws>xw) && (state.T>273.15))
                h = Modelica.Media.Water.IF97_Utilities.h_pT(state.p, state.T, region=1)
            else
                h = 0
            end
            
            return 
        end
        function enthalpyOfGas(
            #= @extends enthalpyOfGas
             =#
            )
            xw::Float64,                  
            xws::Float64,                 
            pd::Float64,                  
            pl::Float64,                  
            pd = Utilities.pd_pTX(state.p, state.T, state.X)
                        pl = state.p-pd
                        xw = state.X[1]/(1-state.X[1])
                        xws = Utilities.xws_pT(state.p, state.T)
                        if ((xw<=xws) || (xws==-1))
                h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, state.T)+xw*Utilities.IF97_new.h_pT(pd, state.T)
            else
                if (state.T<273.16)
                h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, state.T)+xws*Utilities.IF97_new.h_pT(pd, state.T)
            else
                h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, state.T)+xws*Utilities.IF97_new.h_pT(pd, state.T)
            end
            end
            
            return 
        end
        function enthalpyOfCondensingGas(
            #= @extends enthalpyOfCondensingGas
             =#
            )
            xw::Float64,                  
            pd::Float64,                  
            pd = Utilities.pd_pTX(state.p, state.T, state.X)
                        xw = state.X[1]/(1-state.X[1])
                        h = xw*Utilities.IF97_new.h_pT(pd, state.T)
            
            return 
        end
        function enthalpyOfNonCondensingGas(
            #= @extends enthalpyOfNonCondensingGas
             =#
            )
            h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(state.p, state.T)
            
            return 
        end

        "Return specific enthalpy of water (solid + liquid + steam)"
        function enthalpyOfWater(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns h::Float64,            # Specific enthalpy of water
            )
            h = specificEnthalpy(state)-enthalpyOfNonCondensingGas(state)
            
            return h
        end

        "Return enthalpy of liquid and solid water"
        function enthalpyOfWaterNonVapor(
            #=
            @extends Modelica.Icons.Function()
            =#
            state::ThermodynamicState,       # Thermodynamic state record
            # returns h::Float64,            # Specific enthalpy of water
            )
            h = enthalpyOfWater(state)-enthalpyOfWaterVapor(state)
            
            return h
        end
        function pressure(
            #= @extends pressure
             =#
            )
            p = state.p
            
            return 
        end
        function temperature(
            #= @extends temperature
             =#
            )
            T = state.T
            
            return 
        end
        function density(
            #= @extends density
             =#
            )
            d = Utilities.rho_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificEnthalpy(
            #= @extends specificEnthalpy
             =#
            )
            h = Modelica.Media.Air.ReferenceMoistAir.Utilities.h_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificInternalEnergy(
            #= @extends specificInternalEnergy
             =#
            )
            u = Utilities.u_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificEntropy(
            #= @extends specificEntropy
             =#
            )
            s = Utilities.s_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificGibbsEnergy(
            #= @extends specificGibbsEnergy
             =#
            )
            g = Utilities.h_pTX(state.p, state.T, state.X)-state.T*Utilities.s_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificHelmholtzEnergy(
            #= @extends specificHelmholtzEnergy
             =#
            )
            f = Utilities.u_pTX(state.p, state.T, state.X)-state.T*Utilities.s_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificHeatCapacityCp(
            #= @extends specificHeatCapacityCp
             =#
            )
            cp = Utilities.cp_pTX(state.p, state.T, state.X)
            
            return 
        end
        function specificHeatCapacityCv(
            #= @extends specificHeatCapacityCv
             =#
            )
            cv = Utilities.cv_pTX(state.p, state.T, state.X)
            
            return 
        end
        function isentropicExponent(
            #= @extends isentropicExponent
             =#
            )
            gamma = specificHeatCapacityCp(state)/specificHeatCapacityCv(state)
            
            return 
        end
        function isentropicEnthalpy(
            #= @extends isentropicEnthalpy
             =#
            )
            X::Array{MassFraction, 1},       # Complete X-vector
            X = refState.X
                        h_is = specificEnthalpy(setState_psX(p_downstream, specificEntropy(refState), X))
            
            return 
        end
        function velocityOfSound(
            #= @extends velocityOfSound
             =#
            )
            a = sqrt(max(0, gasConstant(state)*state.T*specificHeatCapacityCp(state)/specificHeatCapacityCv(state)))
            
            return 
        end
        function molarMass(
            #= @extends molarMass
             =#
            )
            MM = 1/
                (state.X[1]*steam.MM+state.X[2]*dryair.MM)
            
            return 
        end
        function dynamicViscosity(
            #= @extends dynamicViscosity
             =#
            )
            eta = Utilities.Transport.eta_pTX(state.p, state.T, state.X)
            
            return 
        end
        function thermalConductivity(
            #= @extends thermalConductivity
             =#
            )
            lambda = Utilities.Transport.lambda_pTX(state.p, state.T, state.X)
            
            return 
        end

        "Utility package for moist air"
        module Utilities
            #=
            @extends Modelica.Icons.UtilitiesPackage()
            =#
            const MMX = Array{MoleFraction, 1}=[18.015257E-003,28.01348E-003,31.9988E-003,39.948E-003]
            const Xi_Air = Missing        

            "Compute inverse function"
            module Inverses
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Return temperature as a function of pressure, specific enthalpy and mass fractions"
                function T_phX(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    h::Float64,                      # Specific enthalpy
                    X::Array{Float64, 1},            # Mass fractions
                    # returns T::Float64,            # Temperature
                    )
                    Xfull::Array{MassFraction, 1}=if size(X, 1)==nX; X else cat(1, X, [1-sum(X)]) end,
#=

                    
                    function T_phX_res(
                        #=
                        @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                        =#
                        p::Float64,                      # Pressure
                        h::Float64,                      # Specific enthalpy
                        X::Array{Float64, 1},            # Mass fractions
                        )
                        y = Modelica.Media.Air.ReferenceMoistAir.Utilities.h_pTX(p=p, T=u, X=X)-h
                        
                        return 
                    end
=#
                    T = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
                    
                    return T
                end

                "Return temperature as function of pressure, specific entropy and mass fractions"
                function T_psX(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    s::Float64,                      # Specific entropy
                    X::Array{Float64, 1},            # Mass fractions
                    # returns T::Float64,            # Temperature
                    )
                    Xfull::Array{MassFraction, 1}=if size(X, 1)==nX; X else cat(1, X, [1-sum(X)]) end,
#=

                    
                    function T_psX_res(
                        #=
                        @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                        =#
                        p::Float64,                      # Pressure
                        s::Float64,                      # Specific entropy
                        X::Array{Float64, 1},            # Mass fractions
                        )
                        y = Modelica.Media.Air.ReferenceMoistAir.Utilities.s_pTX(p=p, T=u, X=X)-s
                        
                        return 
                    end
=#
                    T = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
                    
                    return T
                end

                "Return pressure as function of density, temperature and mass fractions"
                function p_dTX(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    X::Array{Float64, 1},            # Mass fractions
                    # returns p::Float64,            # Pressure
                    )
                    Xfull::Array{MassFraction, 1}=if size(X, 1)==nX; X else cat(1, X, [1-sum(X)]) end,
#=

                    
                    function p_dTX_res(
                        #=
                        @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                        =#
                        d::Float64,                      # Density
                        T::Float64,                      # Temperature
                        X::Array{Float64, 1},            # Mass fractions
                        )
                        y = Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX(p=u, T=T, X=X)-d
                        
                        return 
                    end
=#
                    p = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
                    
                    return p
                end
            end

            "Package for transport properties of moist air"
            module Transport
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Coefficients for polynomials used to calculate transport properties"
                mutable struct coef
                    #=
                    @extends Modelica.Icons.Record()
                    =#
                    sigma::Float64                
                    epsilon::Float64              
                    M::Float64                    
                    R::Float64                    
                    w::Array{Float64, 1}          
                    a::Array{Float64, 1}          
                end

                "Dynamic viscosity"
                function eta_pTX(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    X::Array{Float64, 1},            # Mass fractions
                    # returns eta::Float64,          # Dynamic viscosity
                    )
                    ya::Float64,                  
                    yd::Float64,                  
                    yf::Float64,                  
                    va::Float64,                  
                    vd::Float64,                  
                    vf::Float64,                  
                    xw::Float64,                  
                    xws::Float64,                 
                    pd::Float64,                  
                    pl::Float64,                  
                    da::Float64,                  
                    dd::Float64,                  
                    df::Float64,                  
                    Omega::Float64,               
                    Tred::Float64,                
                    etad::Float64,                
                    coef::Modelica.Media.Air.ReferenceMoistAir.Utilities.Transport.coef,
                    xw = X[1]/(1-X[1])
                                        xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                                        pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                                        pl = p-pd
                                        da = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)
                                        if ((xw<=xws) || (xws==-1))
                        if (T<273.16)
                        dd = pd/
                        (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)
                        ya = da/(da+dd)
                        yd = 1-ya
                        Tred = T/coef.epsilon
                        Omega = coef.w[1]+coef.w[2]*Tred+coef.w[3]*Modelica.Math.exp(coef.w[4]*Tred)/(coef.w[5]+Tred)
                        etad = 2.6695E-006*sqrt(T*coef.M)/(coef.sigma^2*Omega)
                        eta = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.eta_dT(da, T)+yd*etad
                    else
                        dd = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.rho_pT(pd, T)
                        ya = da/(da+dd)
                        yd = 1-ya
                        eta = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.eta_dT(da, T)+yd*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.visc_dT(dd, T)
                    end
                    else
                        if (T<273.16)
                        dd = pd/
                        (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)
                        ya = da/(da+dd)
                        yd = 1-ya
                        Tred = T/coef.epsilon
                        Omega = coef.w[1]+coef.w[2]*Tred+coef.w[3]*Modelica.Math.exp(coef.w[4]*Tred)/(coef.w[5]+Tred)
                        etad = 2.6695E-006*sqrt(T*coef.M)/(coef.sigma^2*Omega)
                        eta = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.eta_dT(da, T)+yd*etad
                    else
                        dd = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.rho_pT(pd, T)
                        df = Modelica.Media.Water.IF97_Utilities.rho_pT(p, T)
                        yf = (xw-xws)/df/((1+xws)/(da+dd)+(xw-xws)/df)
                        ya = (1-yf)/(1+dd/da)
                        yd = 1-(ya+yf)
                        eta = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.eta_dT(da, T)+yd*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.visc_dT(dd, T)+yf*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.visc_dT(df, T)
                    end
                    end
                    
                    return eta
                end

                "Thermal conductivity"
                function lambda_pTX(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    X::Array{Float64, 1},            # Mass fractions
                    # returns lambda::Float64,       # Thermal conductivity
                    )
                    ya::Float64,                  
                    yd::Float64,                  
                    yf::Float64,                  
                    va::Float64,                  
                    vd::Float64,                  
                    vf::Float64,                  
                    xw::Float64,                  
                    xws::Float64,                 
                    pd::Float64,                  
                    pl::Float64,                  
                    da::Float64,                  
                    dd::Float64,                  
                    df::Float64,                  
                    Omega::Float64,               
                    Tred::Float64,                
                    cp::Float64,                  
                    Eu::Float64,                  
                    lambdad::Float64,             
                    coef::Modelica.Media.Air.ReferenceMoistAir.Utilities.Transport.coef,
                    xw = X[1]/(1-X[1])
                                        xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                                        pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                                        pl = p-pd
                                        da = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)
                                        if ((xw<=xws) || (xws==-1))
                        if (T<273.16)
                        dd = pd/
                        (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)
                        ya = da/(da+dd)
                        yd = 1-ya
                        Tred = T/coef.epsilon
                        Omega = coef.w[1]+coef.w[2]*Tred+coef.w[3]*Modelica.Math.exp(coef.w[4]*Tred)/(coef.w[5]+Tred)
                        cp = coef.a[1]+coef.a[2]*T+coef.a[3]*T^2+coef.a[4]*T^3+coef.a[5]*T^4
                        Eu = 0.35424*cp+0.1144
                        Eu = 0
                        lambdad = 0.083232*sqrt(T/coef.M)/(coef.sigma^2*Omega)*Eu
                        lambda = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.lambda_dT(da, T)+yd*lambdad
                    else
                        dd = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.rho_pT(pd, T)
                        ya = da/(da+dd)
                        yd = 1-ya
                        lambda = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.lambda_dT(da, T)+yd*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cond_dT(dd, T)
                    end
                    else
                        if (T<273.16)
                        dd = pd/
                        (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)
                        df = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT(p, T)
                        yf = (xw-xws)/df/((1+xws)/(da+dd)+(xw-xws)/df)
                        ya = (1-yf)/(1+dd/da)
                        yd = 1-(ya+yf)
                        Tred = T/coef.epsilon
                        Omega = coef.w[1]+coef.w[2]*Tred+coef.w[3]*Modelica.Math.exp(coef.w[4]*Tred)/(coef.w[5]+Tred)
                        cp = coef.a[1]+coef.a[2]*T+coef.a[3]*T^2+coef.a[4]*T^3+coef.a[5]*T^4
                        Eu = 0.35424*cp+0.1144
                        lambdad = 0.083232*sqrt(T/coef.M)/(coef.sigma^2*Omega)*Eu
                        lambda = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.lambda_dT(da, T)+yd*lambdad+yf*2.21
                    else
                        dd = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.rho_pT(pd, T)
                        df = Modelica.Media.Water.IF97_Utilities.rho_pT(p, T)
                        yf = (xw-xws)/df/((1+xws)/(da+dd)+(xw-xws)/df)
                        ya = (1-yf)/(1+dd/da)
                        yd = 1-(ya+yf)
                        lambda = ya*Modelica.Media.Air.ReferenceAir.Air_Utilities.Transport.lambda_dT(da, T)+yd*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cond_dT(dd, T)+yf*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cond_dT(df, T)
                    end
                    end
                    
                    return lambda
                end
            end

            "Virial and cross-virial coefficients of air and water"
            module VirialCoefficients
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Second molar virial coefficient of dry air"
                function Baa_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns baa::Float64,          # Second virial coefficient
                    )
                    const N::Array{Float64, 1} = Missing,
                    const i::Array{Int64, 1} = Missing,
                    const j::Array{Float64, 1} = Missing,
                    tau::Float64,                 
                    baa = 0
                                        for k in 1:19
                    baa = if (i[k]==1); baa+N[k]*tau^j[k] else baa end
                    end
                    
                                        baa = 1/ReferenceAir.Air_Utilities.Basic.Constants.rhored*baa
                    
                    return baa
                end

                "Second molar cross-virial coefficient"
                function Baw_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns baw::Float64,          # Second cross-virial coefficient
                    )
                    const a::Array{Float64, 1} = Missing,
                    const b::Array{Float64, 1} = Missing,
                    theta::Float64,               
                    baw = 0
                                        theta = T/100
                                        for k in 1:3
                    baw = baw+a[k]*theta^b[k]
                    end
                    
                                        baw = baw*1E-006
                    
                    return baw
                end

                "Second molar virial coefficient of water"
                function Bww_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns bww::Float64,          # Second virial coefficient
                    )
                    const N::Array{Float64, 1} = Missing,
                    const c::Array{Int64, 1} = Missing,
                    const dd::Array{Int64, 1} = Missing,
                    const t::Array{Float64, 1} = Missing,
                    const alpha::Array{Int64, 1} = Missing,
                    const beta::Array{Float64, 1} = Missing,
                    const gamma::Array{Float64, 1} = Missing,
                    const epsilon::Array{Int64, 1} = Missing,
                    const a::Array{Float64, 1} = Missing,
                    const b::Array{Float64, 1} = Missing,
                    const AA::Array{Float64, 1} = Missing,
                    const BB::Array{Float64, 1} = Missing,
                    const CC::Array{Int64, 1} = Missing,
                    const DD::Array{Int64, 1} = Missing,
                    Delta55::Float64,             
                    Delta56::Float64,             
                    theta55::Float64,             
                    theta56::Float64,             
                    psi55::Float64,               
                    psi56::Float64,               
                    Delta55delta::Float64,        
                    Delta56delta::Float64,        
                    Deltab55delta::Float64,       
                    Deltab56delta::Float64,       
                    psi55delta::Float64,          
                    psi56delta::Float64,          
                    tau::Float64,                 
                    bww = 0
                                        theta55 = (1-tau)+AA[55]
                                        theta56 = (1-tau)+AA[56]
                                        psi55 = exp(-CC[55]-DD[55]*(tau-1)^2)
                                        psi56 = exp(-CC[56]-DD[56]*(tau-1)^2)
                                        Delta55 = theta55^2+BB[55]
                                        Delta56 = theta56^2+BB[56]
                                        Delta55delta = -
                        (AA[55]*theta55*2/beta[55]+2*BB[55]*a[55])
                                        Delta56delta = -
                        (AA[56]*theta56*2/beta[56]+2*BB[56]*a[56])
                                        Deltab55delta = b[55]*Delta55^(b[55]-1)*Delta55delta
                                        Deltab56delta = b[56]*Delta56^(b[56]-1)*Delta56delta
                                        psi55delta = 2*CC[55]*psi55
                                        psi56delta = 2*CC[56]*psi56
                                        for k in 1:7
                    bww = if (dd[k]==1); bww+N[k]*tau^t[k] else bww end
                    end
                    
                                        for k in 8:51
                    bww = if (dd[k]==1); bww+N[k]*tau^t[k]*dd[k] else bww end
                    end
                    
                                        for k in 52:54
                    bww = if (dd[k]==1); bww+N[k]*tau^t[k]*exp(-alpha[k]*epsilon[k]^2-beta[k]*(tau-gamma[k])^2) else bww end
                    end
                    
                                        bww = 
                        (bww+N[55]*Delta55^b[55]*psi55+N[56]*Delta56^b[56]*psi56)/Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Constants.rhored*Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Constants.MM
                    
                    return bww
                end

                "Third molar virial coefficient of dry air"
                function Caaa_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns caaa::Float64,         # Third virial coefficient
                    )
                    const N::Array{Float64, 1} = Missing,
                    const i::Array{Int64, 1} = Missing,
                    const j::Array{Float64, 1} = Missing,
                    const l::Array{Int64, 1} = Missing,
                    tau::Float64,                 
                    caaa = 0
                                        for k in 1:10
                    caaa = if (i[k]==2); caaa+2*N[k]*tau^j[k] else caaa end
                    end
                    
                                        for k in 11:19
                    caaa = if (i[k]==2); caaa+2*N[k]*tau^j[k] elseif ((i[k]==1) && (l[k]==1)); caaa-2*N[k]*tau^j[k] else caaa end
                    end
                    
                                        caaa = 1/ReferenceAir.Air_Utilities.Basic.Constants.rhored^2*caaa
                    
                    return caaa
                end

                "Third molar cross-virial coefficient"
                function Caaw_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns caaw::Float64,         # Third cross-virial coefficient
                    )
                    const c::Array{Float64, 1} = Missing,
                    theta::Float64,               
                    caaw = 0
                                        for k in 1:5
                    caaw = caaw+c[k]*theta^(1-k)
                    end
                    
                                        caaw = caaw*1E-012
                    
                    return caaw
                end

                "Third molar cross-virial coefficient"
                function Caww_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns caww::Float64,         # Third cross-virial coefficient
                    )
                    const dd::Array{Float64, 1} = Missing,
                    theta::Float64,               
                    caww = 0
                                        for k in 1:4
                    caww = caww+dd[k]*theta^(1-k)
                    end
                    
                                        caww = -exp(caww)*1E-012
                    
                    return caww
                end

                "Third molar virial coefficient of water"
                function Cwww_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature
                    # returns cwww::Float64,         # Third virial coefficient
                    )
                    const N::Array{Float64, 1} = Missing,
                    const c::Array{Int64, 1} = Missing,
                    const dd::Array{Int64, 1} = Missing,
                    const t::Array{Float64, 1} = Missing,
                    const alpha::Array{Int64, 1} = Missing,
                    const beta::Array{Float64, 1} = Missing,
                    const gamma::Array{Float64, 1} = Missing,
                    const epsilon::Array{Int64, 1} = Missing,
                    const a::Array{Float64, 1} = Missing,
                    const b::Array{Float64, 1} = Missing,
                    const AA::Array{Float64, 1} = Missing,
                    const BB::Array{Float64, 1} = Missing,
                    const CC::Array{Int64, 1} = Missing,
                    const DD::Array{Int64, 1} = Missing,
                    Delta55::Float64,             
                    Delta56::Float64,             
                    theta55::Float64,             
                    theta56::Float64,             
                    psi55::Float64,               
                    psi56::Float64,               
                    Delta55delta::Float64,        
                    Delta56delta::Float64,        
                    Deltab55delta::Float64,       
                    Deltab56delta::Float64,       
                    psi55delta::Float64,          
                    psi56delta::Float64,          
                    Delta55deltadelta::Float64,   
                    Delta56deltadelta::Float64,   
                    Deltab55deltadelta::Float64,  
                    Deltab56deltadelta::Float64,  
                    psi55deltadelta::Float64,     
                    psi56deltadelta::Float64,     
                    tau::Float64,                 
                    cwww = 0
                                        theta55 = (1-tau)+AA[55]
                                        theta56 = (1-tau)+AA[56]
                                        psi55 = exp(-CC[55]-DD[55]*(tau-1)^2)
                                        psi56 = exp(-CC[56]-DD[56]*(tau-1)^2)
                                        Delta55 = theta55^2+BB[55]
                                        Delta56 = theta56^2+BB[56]
                                        Delta55delta = -
                        (AA[55]*theta55*2/beta[55]+2*BB[55]*a[55])
                                        Delta56delta = -
                        (AA[56]*theta56*2/beta[56]+2*BB[56]*a[56])
                                        Deltab55delta = b[55]*Delta55^(b[55]-1)*Delta55delta
                                        Deltab56delta = b[56]*Delta56^(b[56]-1)*Delta56delta
                                        psi55delta = 2*CC[55]*psi55
                                        psi56delta = 2*CC[56]*psi56
                                        Delta55deltadelta = -Delta55delta+AA[55]^2*2/beta[55]^2+AA[55]*theta55*4/beta[55]*(1/(2*beta[55])-1)+4*BB[55]*a[55]*(a[55]-1)
                                        Delta56deltadelta = -Delta56delta+AA[56]^2*2/beta[56]^2+AA[56]*theta56*4/beta[56]*(1/(2*beta[56])-1)+4*BB[56]*a[56]*(a[56]-1)
                                        Deltab55deltadelta = b[55]*
                        (Delta55^(b[55]-1)*Delta55deltadelta+(b[55]-1)*Delta55^(b[55]-2)*Delta55delta^2)
                                        Deltab56deltadelta = b[56]*
                        (Delta56^(b[56]-1)*Delta56deltadelta+(b[56]-1)*Delta56^(b[56]-2)*Delta56delta^2)
                                        psi55deltadelta = (2*CC[55]-1)*2*CC[55]*psi55
                                        psi56deltadelta = (2*CC[56]-1)*2*CC[56]*psi56
                                        cwww = 0
                                        for k in 1:7
                    cwww = if (dd[k]==2); cwww+2*N[k]*tau^t[k] else cwww end
                    end
                    
                                        for k in 8:51
                    cwww = if (dd[k]==2); cwww+2*N[k]*tau^t[k] elseif ((dd[k]==1) && (c[k]==1)); cwww-2*N[k]*tau^t[k] else cwww end
                    end
                    
                                        cwww = cwww+N[55]*
                        (Delta55^b[55]*2*psi55delta+2*Deltab55delta*psi55)+N[56]*
                        (Delta56^b[56]*2*psi56delta+2*Deltab56delta*psi56)
                                        cwww = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Constants.MM^2/Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Constants.rhored^2*cwww*1E-006
                    
                    return cwww
                end
            end

            "Parameters and equations for determining reaction variables (dissociation VDI 4670)"
            module ReactionIndices
                #=
                @extends Modelica.Icons.BasesPackage()
                =#
                const AA = Missing            
                const BB = Missing            
                const CC = Missing            
                const DD = Missing            
                const EE = Missing            
                const p0 = Missing               # Reference pressure

                "Reaction index for formation of H2"
                function U2(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    # returns u::Float64,            # Reaction index for H2
                    )
                    u = AA[2]*moleFraction[1]/sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[2]/T)
                    
                    return u
                end

                "Reaction index for formation of OH"
                function U3(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    # returns u::Float64,            # Reaction index for OH
                    )
                    u = AA[3]*sqrt(moleFraction[1])*sqrt(sqrt(moleFraction[3]))*(p/p0)^(-0.25)*Modelica.Math.exp(BB[3]/T)
                    
                    return u
                end

                "Reaction index for formation of H"
                function U4(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    # returns u::Float64,            # Reaction index for H
                    )
                    u = AA[4]*sqrt(U2(p, T, moleFraction))*(p/p0)^(-0.5)*Modelica.Math.exp(BB[4]/T)
                    
                    return u
                end

                "Reaction index for formation of O"
                function U5(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    # returns u::Float64,            # Reaction index for O
                    )
                    u = AA[5]*sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[5]/T)
                    
                    return u
                end

                "Reaction index for formation of NO"
                function U6(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    # returns u::Float64,            # Reaction index for NO
                    )
                    u = AA[6]*sqrt(moleFraction[2]*moleFraction[3])*Modelica.Math.exp(BB[6]/T)
                    
                    return u
                end

                "Energy index for formation of H2"
                function V2(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns v::Float64,            # Energy index for H2
                    )
                    v = CC[2]+DD[2]/T+EE[2]/T^2
                    
                    return v
                end

                "Energy index for formation of OH"
                function V3(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns v::Float64,            # Energy index for OH
                    )
                    v = CC[3]+DD[3]/T+EE[3]/T^2
                    
                    return v
                end

                "Energy index for formation of H"
                function V4(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns v::Float64,            # Energy index for H
                    )
                    v = CC[4]+DD[4]/T+EE[4]/T^2
                    
                    return v
                end

                "Energy index for formation of O"
                function V5(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns v::Float64,            # Energy index for O
                    )
                    v = CC[5]+DD[5]/T+EE[5]/T^2
                    
                    return v
                end

                "Energy index for formation of NO"
                function V6(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns v::Float64,            # Energy index for NO
                    )
                    v = CC[6]+DD[6]/T+EE[6]/T^2
                    
                    return v
                end

                "Derivative reaction index for formation of H2"
                function U2_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    moleFraction_der::Array{Float64, 1},   # Derivative of mole fractions
                    # returns u_der::Float64,        # Derivative of reaction index for H2
                    )
                    o::Array{Float64, 1},         
                    o[1] = AA[2]*sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[2]/T)
                                        o[2] = -0.5*AA[2]*moleFraction[1]/(moleFraction[3])^1.5*(p/p0)^(-0.5)*Modelica.Math.exp(BB[2]/T)
                                        o[3] = -0.5*AA[2]*moleFraction[1]/sqrt(moleFraction[3])*sqrt(p0)*p^(-1.5)*Modelica.Math.exp(BB[2]/T)
                                        o[4] = -BB[2]*AA[2]*moleFraction[1]/sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[2]/T)/T^2
                                        u_der = o[1]*moleFraction_der[1]+o[2]*moleFraction_der[3]+o[3]*p_der+o[4]*T_der
                    
                    return u_der
                end

                "Derivative of reaction index for formation of OH"
                function U3_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    moleFraction_der::Array{Float64, 1},   # Derivative of mole fractions
                    # returns u_der::Float64,        # Derivative of reaction index for OH
                    )
                    o::Array{Float64, 1},         
                    o[1] = 0.5*AA[3]/sqrt(moleFraction[1])*sqrt(sqrt(moleFraction[3]))*(p/p0)^(-0.25)*Modelica.Math.exp(BB[3]/T)
                                        o[2] = 0.25*AA[3]*sqrt(moleFraction[1])/(moleFraction[3])^0.75*(p/p0)^(-0.25)*Modelica.Math.exp(BB[3]/T)
                                        o[3] = -0.25*AA[3]*sqrt(moleFraction[1])*sqrt(sqrt(moleFraction[3]))*sqrt(sqrt(p0))*p^(-1.25)*Modelica.Math.exp(BB[3]/T)
                                        o[4] = BB[3]*AA[3]*moleFraction[1]/sqrt(sqrt(moleFraction[3]))*(p/p0)^(-0.25)*Modelica.Math.exp(BB[3]/T)/T^2
                                        u_der = o[1]*moleFraction_der[1]+o[2]*moleFraction_der[3]+o[3]*p_der+o[4]*T_der
                    
                    return u_der
                end

                "Derivative of reaction index for formation of H"
                function U4_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    moleFraction_der::Array{Float64, 1},   # Derivative of mole fractions
                    # returns u_der::Float64,        # Derivative of reaction index for H
                    )
                    o::Array{Float64, 1},         
                    o[1] = 0.5*AA[4]/sqrt(U2(p, T, moleFraction))*U2_der(p, T, moleFraction, p_der, T_der, moleFraction_der)*(p/p0)^(-0.5)*Modelica.Math.exp(BB[4]/T)
                                        o[2] = -0.5*AA[4]*sqrt(U2(p, T, moleFraction))*sqrt(p0)*p^(-1.5)*Modelica.Math.exp(BB[4]/T)
                                        o[3] = BB[4]*AA[4]*sqrt(U2(p, T, moleFraction))*(p/p0)^(-0.5)*Modelica.Math.exp(BB[4]/T)/T^2
                                        u_der = o[1]+o[2]*p_der+o[3]*T_der
                    
                    return u_der
                end

                "Derivative of reaction index for formation of O"
                function U5_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    moleFraction_der::Array{Float64, 1},   # Derivative of mole fractions
                    # returns u_der::Float64,        # Derivative of reaction index for O
                    )
                    o::Array{Float64, 1},         
                    o[1] = 0.5*AA[5]/sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[5]/T)
                                        o[2] = -0.5*AA[5]*sqrt(moleFraction[3])*sqrt(p0)*p^(-1.5)*Modelica.Math.exp(BB[5]/T)
                                        o[3] = BB[5]*AA[5]*sqrt(moleFraction[3])*(p/p0)^(-0.5)*Modelica.Math.exp(BB[5]/T)/T^2
                                        u_der = o[1]*moleFraction[3]+o[2]*p_der+o[3]*T_der
                    
                    return u_der
                end

                "Derivative of reaction index for formation of NO"
                function U6_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    moleFraction::Array{Float64, 1},   # Mole fractions
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    moleFraction_der::Array{Float64, 1},   # Derivative of mole fractions
                    # returns u_der::Float64,        # Derivative of reaction index for NO
                    )
                    o::Array{Float64, 1},         
                    o[1] = 0.5*AA[6]/sqrt(moleFraction[2]*moleFraction[3])*moleFraction[3]*Modelica.Math.exp(BB[6]/T)
                                        o[2] = 0.5*AA[6]/sqrt(moleFraction[2]*moleFraction[3])*moleFraction[2]*Modelica.Math.exp(BB[6]/T)
                                        o[3] = BB[6]*AA[6]*sqrt(moleFraction[2]*moleFraction[3])*Modelica.Math.exp(BB[6]/T)/T^2
                                        u_der = o[1]*moleFraction[2]+o[2]*moleFraction[3]+o[3]*T_der
                    
                    return u_der
                end

                "Derivative of energy index for formation of H2"
                function V2_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns v_der::Float64,        # Derivative energy index for H2
                    )
                    v_der = (-DD[2]/T^2-2*EE[2]/T^3)*T_der
                    
                    return v_der
                end

                "Derivative energy index for formation of OH"
                function V3_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns v_der::Float64,        # Derivative energy index for OH
                    )
                    v_der = (-DD[3]/T^2-2*EE[3]/T^3)*T_der
                    
                    return v_der
                end

                "Derivative of energy index for formation of H"
                function V4_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns v_der::Float64,        # Derivative energy index for H
                    )
                    v_der = (-DD[4]/T^2-2*EE[4]/T^3)*T_der
                    
                    return v_der
                end

                "Derivative of energy index for formation of O"
                function V5_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns v_der::Float64,        # Derivative energy index for O
                    )
                    v_der = (-DD[5]/T^2-2*EE[5]/T^3)*T_der
                    
                    return v_der
                end

                "Derivative of energy index for formation of NO"
                function V6_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns v_der::Float64,        # Derivative energy index for NO
                    )
                    v_der = (-DD[6]/T^2-2*EE[6]/T^3)*T_der
                    
                    return v_der
                end
            end

            "Workaround for IF97"
            module IF97_new
                #=
                @extends Modelica.Icons.BasesPackage()
                =#
                const molarMass = MolarMass=0.018015257

                "Gibbs function for region 2: g(p,T)"
                function g2(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature (K)
                    # returns g::Modelica.Media.Common.GibbsDerivs,   # Dimensionless Gibbs function and dervatives w.r.t. pi and tau
                    )
                    tau2::Float64,                   # Dimensionless temperature
                    o::Array{Float64, 1},            # Vector of auxiliary variables
                    g.p = p
                                        g.T = T
                                        g.R = Modelica.Media.Water.IF97_Utilities.BaseIF97.data.RH2O
                                        g.pi = p/Modelica.Media.Water.IF97_Utilities.BaseIF97.data.PSTAR2
                                        g.tau = Modelica.Media.Water.IF97_Utilities.BaseIF97.data.TSTAR2/T
                                        tau2 = -0.5+g.tau
                                        o[1] = tau2*tau2
                                        o[2] = o[1]*tau2
                                        o[3] = -0.050325278727930*o[2]
                                        o[4] = -0.057581259083432+o[3]
                                        o[5] = o[4]*tau2
                                        o[6] = -0.045996013696365+o[5]
                                        o[7] = o[6]*tau2
                                        o[8] = -0.0178348622923580+o[7]
                                        o[9] = o[8]*tau2
                                        o[10] = o[1]*o[1]
                                        o[11] = o[10]*o[10]
                                        o[12] = o[11]*o[11]
                                        o[13] = o[10]*o[11]*o[12]*tau2
                                        o[14] = o[1]*o[10]*tau2
                                        o[15] = o[10]*o[11]*tau2
                                        o[16] = o[1]*o[12]*tau2
                                        o[17] = o[1]*o[11]*tau2
                                        o[18] = o[1]*o[10]*o[11]
                                        o[19] = o[10]*o[11]*o[12]
                                        o[20] = o[1]*o[10]
                                        o[21] = g.pi*g.pi
                                        o[22] = o[21]*o[21]
                                        o[23] = o[21]*o[22]
                                        o[24] = o[10]*o[12]*tau2
                                        o[25] = o[12]*o[12]
                                        o[26] = o[11]*o[12]*o[25]*tau2
                                        o[27] = o[10]*o[12]
                                        o[28] = o[1]*o[10]*o[11]*tau2
                                        o[29] = o[10]*o[12]*o[25]*tau2
                                        o[30] = o[1]*o[10]*o[25]*tau2
                                        o[31] = o[1]*o[11]*o[12]
                                        o[32] = o[1]*o[12]
                                        o[33] = g.tau*g.tau
                                        o[34] = o[33]*o[33]
                                        o[35] = -0.000053349095828174*o[13]
                                        o[36] = -0.087594591301146+o[35]
                                        o[37] = o[2]*o[36]
                                        o[38] = -0.0078785554486710+o[37]
                                        o[39] = o[1]*o[38]
                                        o[40] = -0.00037897975032630+o[39]
                                        o[41] = o[40]*tau2
                                        o[42] = -0.000066065283340406+o[41]
                                        o[43] = o[42]*tau2
                                        o[44] = 5.7870447262208e-6*tau2
                                        o[45] = -0.301951672367580*o[2]
                                        o[46] = -0.172743777250296+o[45]
                                        o[47] = o[46]*tau2
                                        o[48] = -0.091992027392730+o[47]
                                        o[49] = o[48]*tau2
                                        o[50] = o[1]*o[11]
                                        o[51] = o[10]*o[11]
                                        o[52] = o[11]*o[12]*o[25]
                                        o[53] = o[10]*o[12]*o[25]
                                        o[54] = o[1]*o[10]*o[25]
                                        o[55] = o[11]*o[12]*tau2
                                        g.g = g.pi*
                        (-0.00177317424732130+o[9]+g.pi*
                        (tau2*
                        (-0.000033032641670203+
                        (-0.000189489875163150+o[1]*
                        (-0.0039392777243355+
                        (-0.043797295650573-0.0000266745479140870*o[13])*o[2]))*tau2)+g.pi*
                        (2.04817376923090e-8+
                        (4.3870667284435e-7+o[1]*
                        (-0.000032277677238570+
                        (-0.00150339245421480-0.040668253562649*o[13])*o[2]))*tau2+g.pi*
                        (g.pi*
                        (2.29220763376610e-6*o[14]+g.pi*
                        (
                        (-1.67147664510610e-11+o[15]*
                        (-0.00211714723213550-23.8957419341040*o[16]))*o[2]+g.pi*
                        (-5.9059564324270e-18+o[17]*
                        (-1.26218088991010e-6-0.038946842435739*o[18])+g.pi*
                        (o[11]*
                        (1.12562113604590e-11-8.2311340897998*o[19])+g.pi*
                        (1.98097128020880e-8*o[15]+g.pi*
                        (o[10]*
                        (1.04069652101740e-19+
                        (-1.02347470959290e-13-1.00181793795110e-9*o[10])*o[20])+o[23]*
                        (o[13]*
                        (-8.0882908646985e-11+0.106930318794090*o[24])+o[21]*
                        (-0.33662250574171*o[26]+o[21]*
                        (o[27]*
                        (8.9185845355421e-25+
                        (3.06293168762320e-13-4.2002467698208e-6*o[15])*o[28])+g.pi*
                        (-5.9056029685639e-26*o[24]+g.pi*
                        (3.7826947613457e-6*o[29]+g.pi*
                        (-1.27686089346810e-15*o[30]+o[31]*
                        (7.3087610595061e-29+o[18]*
                        (5.5414715350778e-17-9.4369707241210e-7*o[32]))*g.pi))))))))))))+tau2*
                        (-7.8847309559367e-10+
                        (1.27907178522850e-8+4.8225372718507e-7*tau2)*tau2)))))+
                        (-0.0056087911830200+g.tau*
                        (0.071452738814550+g.tau*
                        (-0.40710498239280+g.tau*
                        (1.42408197144400+g.tau*
                        (-4.3839511194500+g.tau*
                        (-9.6927686002170+g.tau*
                        (10.0866556801800+
                        (-0.284086326077200+0.0212684635330700*g.tau)*g.tau)+Modelica.Math.log(g.pi)))))))/(o[34]*g.tau)
                                        g.gpi = 
                        (1.00000000000000+g.pi*
                        (-0.00177317424732130+o[9]+g.pi*
                        (o[43]+g.pi*
                        (6.1445213076927e-8+
                        (1.31612001853305e-6+o[1]*
                        (-0.000096833031715710+
                        (-0.0045101773626444-0.122004760687947*o[13])*o[2]))*tau2+g.pi*
                        (g.pi*
                        (0.0000114610381688305*o[14]+g.pi*
                        (
                        (-1.00288598706366e-10+o[15]*
                        (-0.0127028833928130-143.374451604624*o[16]))*o[2]+g.pi*
                        (-4.1341695026989e-17+o[17]*
                        (-8.8352662293707e-6-0.272627897050173*o[18])+g.pi*
                        (o[11]*
                        (9.0049690883672e-11-65.849072718398*o[19])+g.pi*
                        (1.78287415218792e-7*o[15]+g.pi*
                        (o[10]*
                        (1.04069652101740e-18+
                        (-1.02347470959290e-12-1.00181793795110e-8*o[10])*o[20])+o[23]*
                        (o[13]*
                        (-1.29412653835176e-9+1.71088510070544*o[24])+o[21]*
                        (-6.0592051033508*o[26]+o[21]*
                        (o[27]*
                        (1.78371690710842e-23+
                        (6.1258633752464e-12-0.000084004935396416*o[15])*o[28])+g.pi*
                        (-1.24017662339842e-24*o[24]+g.pi*
                        (0.000083219284749605*o[29]+g.pi*
                        (-2.93678005497663e-14*o[30]+o[31]*
                        (1.75410265428146e-27+o[18]*
                        (1.32995316841867e-15-0.0000226487297378904*o[32]))*g.pi))))))))))))+tau2*
                        (-3.15389238237468e-9+
                        (5.1162871409140e-8+1.92901490874028e-6*tau2)*tau2))))))/g.pi
                                        g.gpipi = 
                        (-1.00000000000000+o[21]*
                        (o[43]+g.pi*
                        (1.22890426153854e-7+
                        (2.63224003706610e-6+o[1]*
                        (-0.000193666063431420+
                        (-0.0090203547252888-0.244009521375894*o[13])*o[2]))*tau2+g.pi*
                        (g.pi*
                        (0.000045844152675322*o[14]+g.pi*
                        (
                        (-5.0144299353183e-10+o[15]*
                        (-0.063514416964065-716.87225802312*o[16]))*o[2]+g.pi*
                        (-2.48050170161934e-16+o[17]*
                        (-0.000053011597376224-1.63576738230104*o[18])+g.pi*
                        (o[11]*
                        (6.3034783618570e-10-460.94350902879*o[19])+g.pi*
                        (1.42629932175034e-6*o[15]+g.pi*
                        (o[10]*
                        (9.3662686891566e-18+
                        (-9.2112723863361e-12-9.0163614415599e-8*o[10])*o[20])+o[23]*
                        (o[13]*
                        (-1.94118980752764e-8+25.6632765105816*o[24])+o[21]*
                        (-103.006486756963*o[26]+o[21]*
                        (o[27]*
                        (3.3890621235060e-22+
                        (1.16391404129682e-10-0.00159609377253190*o[15])*o[28])+g.pi*
                        (-2.48035324679684e-23*o[24]+g.pi*
                        (0.00174760497974171*o[29]+g.pi*
                        (-6.4609161209486e-13*o[30]+o[31]*
                        (4.0344361048474e-26+o[18]*
                        (3.05889228736295e-14-0.00052092078397148*o[32]))*g.pi))))))))))))+tau2*
                        (-9.4616771471240e-9+(1.53488614227420e-7+o[44])*tau2)))))/o[21]
                                        g.gtau = 
                        (0.0280439559151000+g.tau*
                        (-0.285810955258200+g.tau*
                        (1.22131494717840+g.tau*
                        (-2.84816394288800+g.tau*
                        (4.3839511194500+o[33]*
                        (10.0866556801800+
                        (-0.56817265215440+0.063805390599210*g.tau)*g.tau))))))/(o[33]*o[34])+g.pi*
                        (-0.0178348622923580+o[49]+g.pi*
                        (-0.000033032641670203+
                        (-0.00037897975032630+o[1]*
                        (-0.0157571108973420+
                        (-0.306581069554011-0.00096028372490713*o[13])*o[2]))*tau2+g.pi*
                        (4.3870667284435e-7+o[1]*
                        (-0.000096833031715710+
                        (-0.0090203547252888-1.42338887469272*o[13])*o[2])+g.pi*
                        (-7.8847309559367e-10+g.pi*
                        (0.0000160454534363627*o[20]+g.pi*
                        (o[1]*
                        (-5.0144299353183e-11+o[15]*
                        (-0.033874355714168-836.35096769364*o[16]))+g.pi*
                        (
                        (-0.0000138839897890111-0.97367106089347*o[18])*o[50]+g.pi*
                        (o[14]*
                        (9.0049690883672e-11-296.320827232793*o[19])+g.pi*
                        (2.57526266427144e-7*o[51]+g.pi*
                        (o[2]*
                        (4.1627860840696e-19+
                        (-1.02347470959290e-12-1.40254511313154e-8*o[10])*o[20])+o[23]*
                        (o[19]*
                        (-2.34560435076256e-9+5.3465159397045*o[24])+o[21]*
                        (-19.1874828272775*o[52]+o[21]*
                        (o[16]*
                        (1.78371690710842e-23+
                        (1.07202609066812e-11-0.000201611844951398*o[15])*o[28])+g.pi*
                        (-1.24017662339842e-24*o[27]+g.pi*
                        (0.000200482822351322*o[53]+g.pi*
                        (-4.9797574845256e-14*o[54]+
                        (1.90027787547159e-27+o[18]*
                        (2.21658861403112e-15-0.000054734430199902*o[32]))*o[55]*g.pi))))))))))))+
                        (2.55814357045700e-8+1.44676118155521e-6*tau2)*tau2))))
                                        g.gtautau = 
                        (-0.168263735490600+g.tau*
                        (1.42905477629100+g.tau*
                        (-4.8852597887136+g.tau*
                        (8.5444918286640+g.tau*
                        (-8.7679022389000+o[33]*
                        (-0.56817265215440+0.127610781198420*g.tau)*g.tau)))))/(o[33]*o[34]*g.tau)+g.pi*
                        (-0.091992027392730+
                        (-0.34548755450059-1.50975836183790*o[2])*tau2+g.pi*
                        (-0.00037897975032630+o[1]*
                        (-0.047271332692026+
                        (-1.83948641732407-0.033609930371750*o[13])*o[2])+g.pi*
                        (
                        (-0.000193666063431420+
                        (-0.045101773626444-48.395221739552*o[13])*o[2])*tau2+g.pi*
                        (2.55814357045700e-8+2.89352236311042e-6*tau2+g.pi*
                        (0.000096272720618176*o[10]*tau2+g.pi*
                        (
                        (-1.00288598706366e-10+o[15]*
                        (-0.50811533571252-28435.9329015838*o[16]))*tau2+g.pi*
                        (o[11]*
                        (-0.000138839897890111-23.3681054614434*o[18])*tau2+g.pi*
                        (
                        (6.3034783618570e-10-10371.2289531477*o[19])*o[20]+g.pi*
                        (3.09031519712573e-6*o[17]+g.pi*
                        (o[1]*
                        (1.24883582522088e-18+
                        (-9.2112723863361e-12-1.82330864707100e-7*o[10])*o[20])+o[23]*
                        (o[1]*o[11]*o[12]*
                        (-6.5676921821352e-8+261.979281045521*o[24])*tau2+o[21]*
                        (-1074.49903832754*o[1]*o[10]*o[12]*o[25]*tau2+o[21]*
                        (
                        (3.3890621235060e-22+
                        (3.6448887082716e-10-0.0094757567127157*o[15])*o[28])*o[32]+g.pi*
                        (-2.48035324679684e-23*o[16]+g.pi*
                        (0.0104251067622687*o[1]*o[12]*o[25]*tau2+g.pi*
                        (o[11]*o[12]*
                        (4.7506946886790e-26+o[18]*
                        (8.6446955947214e-14-0.00311986252139440*o[32]))*g.pi-1.89230784411972e-12*o[10]*o[25]*tau2))))))))))))))))
                                        g.gtaupi = -0.0178348622923580+o[49]+g.pi*
                        (-0.000066065283340406+
                        (-0.00075795950065260+o[1]*
                        (-0.0315142217946840+
                        (-0.61316213910802-0.00192056744981426*o[13])*o[2]))*tau2+g.pi*
                        (1.31612001853305e-6+o[1]*
                        (-0.000290499095147130+
                        (-0.0270610641758664-4.2701666240781*o[13])*o[2])+g.pi*
                        (-3.15389238237468e-9+g.pi*
                        (0.000080227267181813*o[20]+g.pi*
                        (o[1]*
                        (-3.00865796119098e-10+o[15]*
                        (-0.203246134285008-5018.1058061618*o[16]))+g.pi*
                        (
                        (-0.000097187928523078-6.8156974262543*o[18])*o[50]+g.pi*
                        (o[14]*
                        (7.2039752706938e-10-2370.56661786234*o[19])+g.pi*
                        (2.31773639784430e-6*o[51]+g.pi*
                        (o[2]*
                        (4.1627860840696e-18+
                        (-1.02347470959290e-11-1.40254511313154e-7*o[10])*o[20])+o[23]*
                        (o[19]*
                        (-3.7529669612201e-8+85.544255035272*o[24])+o[21]*
                        (-345.37469089099*o[52]+o[21]*
                        (o[16]*
                        (3.5674338142168e-22+
                        (2.14405218133624e-10-0.0040322368990280*o[15])*o[28])+g.pi*
                        (-2.60437090913668e-23*o[27]+g.pi*
                        (0.0044106220917291*o[53]+g.pi*
                        (-1.14534422144089e-12*o[54]+
                        (4.5606669011318e-26+o[18]*
                        (5.3198126736747e-14-0.00131362632479764*o[32]))*o[55]*g.pi))))))))))))+(1.02325742818280e-7+o[44])*tau2)))
                    
                    return g
                end

                "Specific enthalpy as function or pressure and temperature"
                function h_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    region::Int64,                   # If 0, region is unknown, otherwise known and this input
                    # returns h::Float64,            # Specific enthalpy
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        h = g.R*T*g.tau*g.gtau
                    
                    return h
                end

                "Temperature as function of pressure and temperature"
                function s_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    region::Int64,                   # If 0, region is unknown, otherwise known and this input
                    # returns s::Float64,            # Specific entropy
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        s = g.R*(g.tau*g.gtau-g.g)
                    
                    return s
                end

                "Specific heat capacity at constant pressure as function of pressure and temperature"
                function cp_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    region::Int64,                   # If 0, region is unknown, otherwise known and this input
                    # returns cp::Float64,           # Specific heat capacity
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        cp = -g.R*g.tau*g.tau*g.gtautau
                    
                    return cp
                end

                "Specific heat capacity at constant volume as function of pressure and temperature"
                function cv_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    region::Int64,                   # If 0, region is unknown, otherwise known and this input
                    # returns cv::Float64,           # Specific heat capacity
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        cv = g.R*
                        (-g.tau*g.tau*g.gtautau+
                        ((g.gpi-g.tau*g.gtaupi)*(g.gpi-g.tau*g.gtaupi)/g.gpipi))
                    
                    return cv
                end

                "Density as function or pressure and temperature"
                function rho_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns rho::Float64,          # Density
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        rho = p/(g.R*T*g.pi*g.gpi)
                    
                    return rho
                end

                "Derivative function of rho_pT"
                function rho_pT_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    # returns rho_der::Float64,      # Derivative of density
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    d::Float64,                   
                    vp::Float64,                  
                    vt::Float64,                  
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        vt = g.R/p*(g.pi*g.gpi-g.tau*g.pi*g.gtaupi)
                                        vp = g.R*T/(p*p)*g.pi*g.pi*g.gpipi
                                        d = p/(g.R*T*g.pi*g.gpi)
                                        rho_der = (-d^2*vp)*p_der+(-d^2*vt)*T_der
                    
                    return rho_der
                end

                "Dynamic viscosity eta(d,T), industrial formulation"
                function visc_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature (K)
                    # returns eta::Float64,          # Dynamic viscosity
                    )
                    const n0::Float64 = Missing,     # Viscosity coefficient
                    const n1::Float64 = Missing,     # Viscosity coefficient
                    const n2::Float64 = Missing,     # Viscosity coefficient
                    const n3::Float64 = Missing,     # Viscosity coefficient
                    const nn::Array{Float64, 1} = Missing,   # Viscosity coefficients
                    const rhostar::Float64 = Missing,   # Scaling density
                    const etastar::Float64 = Missing,   # Scaling viscosity
                    const tstar::Float64 = Missing,   # Scaling temperature
                    i::Int64,                        # Auxiliary variable
                    j::Int64,                        # Auxiliary variable
                    delta::Float64,                  # Dimensionless density
                    deltam1::Float64,                # Dimensionless density
                    tau::Float64,                    # Dimensionless temperature
                    taum1::Float64,                  # Dimensionless temperature
                    Psi0::Float64,                   # Auxiliary variable
                    Psi1::Float64,                   # Auxiliary variable
                    tfun::Float64,                   # Auxiliary variable
                    rhofun::Float64,                 # Auxiliary variable
                    delta = d/rhostar
                                        deltam1 = delta-1.0
                                        tau = tstar/T
                                        taum1 = tau-1.0
                                        Psi0 = 1/(n0+(n1+(n2+n3*tau)*tau)*tau)/(tau^0.5)
                                        Psi1 = 0.0
                                        tfun = 1.0
                                        for i in 1:6
                    if (i!=1)
                        tfun = tfun*taum1
                    end
                    rhofun = 1.0
                    for j in 0:6
                    if (j!=0)
                        rhofun = rhofun*deltam1
                    end
                    Psi1 = Psi1+nn[i+j*6]*tfun*rhofun
                    end
                    
                    end
                    
                                        eta = etastar*Psi0*Modelica.Math.exp(delta*Psi1)
                    
                    return eta
                end

                "Thermal conductivity"
                function cond_dT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    d::Float64,                      # Density
                    T::Float64,                      # Temperature (K)
                    # returns lambda::Float64,       # Thermal conductivity
                    )
                    region::Int64,                   # IF97 region, valid values:1,2,3, and 5
                    const n0::Float64 = Missing,     # Conductivity coefficient
                    const n1::Float64 = Missing,     # Conductivity coefficient
                    const n2::Float64 = Missing,     # Conductivity coefficient
                    const n3::Float64 = Missing,     # Conductivity coefficient
                    const nn::Array{Float64, 1} = Missing,   # Conductivity coefficient
                    const lamstar::Float64 = Missing,   # Scaling conductivity
                    const rhostar::Float64 = Missing,   # Scaling density
                    const tstar::Float64 = Missing,   # Scaling temperature
                    const pstar::Float64 = Missing,   # Scaling pressure
                    const etastar::Float64 = Missing,   # Scaling viscosity
                    i::Int64,                        # Auxiliary variable
                    j::Int64,                        # Auxiliary variable
                    delta::Float64,                  # Dimensionless density
                    tau::Float64,                    # Dimensionless temperature
                    deltam1::Float64,                # Dimensionless density
                    taum1::Float64,                  # Dimensionless temperature
                    Lam0::Float64,                   # Part of thermal conductivity
                    Lam1::Float64,                   # Part of thermal conductivity
                    Lam2::Float64,                   # Part of thermal conductivity
                    tfun::Float64,                   # Auxiliary variable
                    rhofun::Float64,                 # Auxiliary variable
                    dpitau::Float64,                 # Auxiliary variable
                    ddelpi::Float64,                 # Auxiliary variable
                    d2::Float64,                     # Auxiliary variable
                    g::Modelica.Media.Common.GibbsDerivs,   # Dimensionless Gibbs function and derivatives w.r.t. pi and tau
                    f::Modelica.Media.Common.HelmholtzDerivs,   # Dimensionless Helmholtz function and derivatives w.r.t. delta and tau
                    Tc::Float64,                     # Celsius temperature for region check
                    Chi::Float64,                    # Symmetrized compressibility
                    const rhostar2::Float64 = Missing,   # Reference density
                    const Tstar2::Float64 = Missing,   # Reference temperature
                    const lambdastar::Float64 = Missing,   # Reference thermal conductivity
                    TREL::Float64,                   # Relative temperature
                    rhoREL::Float64,                 # Relative density
                    lambdaREL::Float64,              # Relative thermal conductivity
                    deltaTREL::Float64,              # Relative temperature increment
                    const C::Array{Float64, 1} = Missing,
                    const dpar::Array{Float64, 1} = Missing,
                    const b::Array{Float64, 1} = Missing,
                    const B::Array{Float64, 1} = Missing,
                    const a::Array{Float64, 1} = Missing,
                    Q::Float64,                   
                    S::Float64,                   
                    lambdaREL2::Float64,             # Function, part of the interpolating equation of the thermal conductivity
                    lambdaREL1::Float64,             # Function, part of the interpolating equation of the thermal conductivity
                    lambdaREL0::Float64,             # Function, part of the interpolating equation of the thermal conductivity
                    deltaTREL = abs(TREL-1)+C[4]
                                        Q = 2+C[5]/deltaTREL^(3/5)
                                        if TREL>=1
                        S = 1/deltaTREL
                    else
                        S = C[6]/deltaTREL^(3/5)
                    end
                                        lambdaREL2 = (dpar[1]/TREL^10+dpar[2])*rhoREL^(9/5)*Modelica.Math.exp(C[1]*(1-rhoREL^(14/5)))+dpar[3]*S*rhoREL^Q*Modelica.Math.exp((Q/(1+Q))*(1-rhoREL^(1+Q)))+dpar[4]*Modelica.Math.exp(C[2]*TREL^(3/2)+C[3]/rhoREL^5)
                                        lambdaREL1 = b[1]+b[2]*rhoREL+b[3]*Modelica.Math.exp(B[1]*(rhoREL+B[2])^2)
                                        lambdaREL0 = TREL^(1/2)*sum(a[i]*TREL^(i-1) for i in 1:4)
                                        lambdaREL = lambdaREL0+lambdaREL1+lambdaREL2
                                        lambda = lambdaREL*lambdastar
                    
                    return lambda
                end

                "Derivative function of h_pT"
                function h_pT_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    region::Int64,                   # If 0, region is unknown, otherwise known and this input
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    # returns h_der::Float64,        # Derivative of specific enthalpy
                    )
                    g::Modelica.Media.Common.GibbsDerivs,
                    rho::Float64,                 
                    vt::Float64,                  
                    g = Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.g2(p, T)
                                        vt = g.R/p*(g.pi*g.gpi-g.tau*g.pi*g.gtaupi)
                                        rho = max(p/(g.R*T*g.pi*g.gpi), 1e-9)
                                        h_der = (1/rho-T*vt)*p_der-g.R*g.tau*g.tau*g.gtautau*T_der
                    
                    return h_der
                end
            end

            "Utility functions from IAPWS95 required for air temperatures below 273.15 K"
            module Water95_Utilities
                #=
                @extends Modelica.Icons.BasesPackage()
                =#
                const Constants = Common.FundamentalConstants(
                    R_bar=8.314371, 
                    R=461.51805, 
                    MM=18.015268E-003, 
                    rhored=322, 
                    Tred=647.096, 
                    pred=22064000, 
                    h_off=0.0, 
                    s_off=0.0)
                const waterConstants = Modelica.Media.Interfaces.Types.TwoPhase.FluidConstants(
                    chemicalFormula="H2O", 
                    structureFormula="H2O", 
                    casRegistryNumber="7732-18-5", 
                    iupacName="oxidane", 
                    molarMass=0.018015268, 
                    criticalTemperature=647.096, 
                    criticalPressure=22064.0e3, 
                    criticalMolarVolume=1/322.0*0.018015268, 
                    triplePointTemperature=273.16, 
                    triplePointPressure=611.657, 
                    normalBoilingPoint=373.124, 
                    meltingPoint=273.15, 
                    acentricFactor=0.344, 
                    dipoleMoment=1.8, 
                    hasCriticalData=true, 
                    hasIdealGasHeatCapacity=false, 
                    hasDipoleMoment=true, 
                    hasFundamentalEquation=true, 
                    hasLiquidHeatCapacity=true, 
                    hasSolidHeatCapacity=false, 
                    hasAccurateViscosityData=true, 
                    hasAccurateConductivityData=true, 
                    hasVapourPressureCurve=false, 
                    hasAcentricFactor=true, 
                    HCRIT0=0.0, 
                    SCRIT0=0.0, 
                    deltah=0.0, 
                    deltas=0.0)

                "Saturation pressure"
                function psat(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    # returns p::Float64,            # Pressure
                    )
                    theta_s::Float64,             
                    A::Float64,                   
                    B::Float64,                   
                    C::Float64,                   
                    const n::Array{Float64, 1} = Missing,
                    theta_s = min(T, 553)+n[9]/(min(T, 553)-n[10])
                                        A = theta_s^2+n[1]*theta_s+n[2]
                                        B = n[3]*theta_s^2+n[4]*theta_s+n[5]
                                        C = n[6]*theta_s^2+n[7]*theta_s+n[8]
                                        p = (2*C/(-B+sqrt(B^2-4*A*C)))^4*1E+006
                    
                    return p
                end

                "Saturation temperature"
                function Tsat(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    # returns T::Float64,            # Temperature
                    )
                    beta::Float64,                
                    D::Float64,                   
                    E::Float64,                   
                    F::Float64,                   
                    G::Float64,                   
                    const n::Array{Float64, 1} = Missing,
                    beta = (p*1E-006)^0.25
                                        E = beta^2+n[3]*beta+n[6]
                                        F = n[1]*beta^2+n[4]*beta+n[7]
                                        G = n[2]*beta^2+n[5]*beta+n[8]
                                        D = 2*G/(-F-sqrt(F^2-4*E*G))
                                        T = 
                        (n[10]+D-sqrt((n[10]+D)^2-4*(n[9]+n[10]*D)))/2
                    
                    return T
                end

                "Saturation pressure"
                function psat_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    T::Float64,                      # Temperature
                    T_der::Float64,                  # Derivative of temperature
                    # returns p_der::Float64,        # Derivative of pressure w.r.t. temperature
                    )
                    theta_s::Float64,             
                    A::Float64,                   
                    B::Float64,                   
                    C::Float64,                   
                    const n::Array{Float64, 1} = Missing,
                    theta_s = min(T, 553)+n[9]/(min(T, 553)-n[10])
                                        theta_s_der = (1-n[9]/(min(T, 553)-n[10])^2)
                                        A = theta_s^2+n[1]*theta_s+n[2]
                                        B = n[3]*theta_s^2+n[4]*theta_s+n[5]
                                        C = n[6]*theta_s^2+n[7]*theta_s+n[8]
                                        A_der = 2*theta_s*theta_s_der+n[1]*theta_s_der
                                        B_der = 2*n[3]*theta_s*theta_s_der+n[4]*theta_s_der
                                        C_der = 2*n[6]*theta_s*theta_s_der+n[7]*theta_s_der
                                        o_der[1] = 2*B*B_der-4*(A*C_der+A_der*C)
                                        o_der[2] = -B_der+0.5/sqrt(B^2-4*A*C)*o_der[1]
                                        o_der[3] = 
                        ((2*C_der*(-B+sqrt(B^2-4*A*C)))-2*C*o_der[2])/(-B+sqrt(B^2-4*A*C))^2
                                        p_der = 4*((2*C/(-B+sqrt(B^2-4*A*C))))^3*o_der[3]*1E+006*T_der
                    
                    return p_der
                end

                "Derivative of saturation temperature"
                function Tsat_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    p_der::Float64,                  # Pressure derivative
                    # returns T_der::Float64,        # Temperature derivative
                    )
                    beta::Float64,                
                    D::Float64,                   
                    E::Float64,                   
                    F::Float64,                   
                    G::Float64,                   
                    o_der::Array{Float64, 1},     
                    const n::Array{Float64, 1} = Missing,
                    beta = (p*1E-006)^0.25
                                        beta_der = 0.25/(p*1E-006)^0.75*1E-006
                                        E = beta^2+n[3]*beta+n[6]
                                        E_der = 2*beta*beta_der+n[3]
                                        F = n[1]*beta^2+n[4]*beta+n[7]
                                        F_der = 2*n[1]*beta*beta_der+n[4]
                                        G = n[2]*beta^2+n[5]*beta+n[8]
                                        G_der = 2*n[2]*beta*beta_der+n[5]
                                        D = 2*G/(-F-sqrt(F^2-4*E*G))
                                        o_der[1] = 2*F*F_der-4*(E*G_der+E_der*G)
                                        o_der[2] = -F_der-0.5/sqrt(F^2-4*E*G)*o_der[1]
                                        D_der = 
                        ((2*G_der*(-F-sqrt(F^2-4*E*G)))-2*G*o_der[2])/(-F-sqrt(F^2-4*E*G))^2
                                        T_der = 
                        (D_der-0.5/sqrt((n[10]+D)^2-4*(n[9]+n[10]*D))*(2*D*D_der-2*n[10]*D_der))/2*p_der
                    
                    return T_der
                end
            end

            "Utility functions from IAPWS09 required for ice in air"
            module Ice09_Utilities
                #=
                @extends Modelica.Icons.BasesPackage()
                =#

                "Fundamental equation of state"
                module Basic
                    #=
                    @extends Modelica.Icons.BasesPackage()
                    =#

                    
                    mutable struct IceConstants
                        #=
                        @extends Common.FundamentalConstants()
                        =#
                        p0::Float64                   
                    end

                    "Complex number with function"
                    mutable struct MyComplex
                        re::Float64                      # Real part of complex number
                        im::Float64                      # Imaginary part of complex number
                    end
                    const Constants = IceConstants(
                        R_bar=8.314472, 
                        R=461.52364, 
                        MM=18.015268E-003, 
                        rhored=1.0, 
                        Tred=273.16, 
                        pred=611.657, 
                        p0=101325, 
                        h_off=0.0, 
                        s_off=0.0)

                    "Sublimation pressure"
                    function psub(
                        #=
                        @extends Modelica.Icons.Function()
                        =#
                        T::Float64,                      # Temperature
                        # returns p_sub::Float64,        # Pressure
                        )
                        const a::Array{Float64, 1} = Missing,
                        const b::Array{Float64, 1} = Missing,
                        theta::Float64,               
                        sum::Float64,                 
                        theta = T/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.Tred
                                                for k in 1:3
                        sum = sum+a[k]*theta^b[k]
                        end
                        
                                                p_sub = exp(sum/theta)*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.pred
                        
                        return p_sub
                    end

                    "Sublimation temperature"
                    function Tsub(
                        #=
                        @extends Modelica.Icons.Function()
                        =#
                        p::Float64,                      # Pressure
                        # returns T_sub::Float64,        # Temperature
                        )
#=

                        
                        function Tsub_res(
                            #=
                            @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                            =#
                            p::Float64,                      # Pressure
                            )
                            y = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub(u)-p
                            
                            return 
                        end
=#
                        T_sub = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
                        
                        return T_sub
                    end

                    "Derivative of sublimation pressure"
                    function psub_der(
                        #=
                        @extends Modelica.Icons.Function()
                        =#
                        T::Float64,                      # Temperature
                        T_der::Float64,                  # Derivative of temperature
                        # returns p_sub_der::Float64,    # Derivative of pressure
                        )
                        const a::Array{Float64, 1} = Missing,
                        const b::Array{Float64, 1} = Missing,
                        theta::Float64,               
                        theta_der::Float64,           
                        sum::Float64,                 
                        sum_der::Float64,             
                        theta = T/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.Tred
                                                theta_der = 1/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.Tred
                                                for k in 1:3
                        sum = sum+a[k]*theta^b[k]
                        sum_der = sum_der+a[k]*b[k]*theta^(b[k]-1)*theta_der
                        end
                        
                                                p_sub_der = (sum_der*theta-sum*theta_der)/theta^2*exp(sum/theta)*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.pred*T_der
                        
                        return p_sub_der
                    end
                end

                "Intermediate property record for water (p and T preferred states)"
                function ice09BaseProp_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns aux::Common.AuxiliaryProperties,   # Auxiliary record
                    )
                    g::Common.GibbsDerivs2,          # Gibbs function and dervatives w.r.t. p and T
                    aux.p = p
                                        aux.T = T
                                        aux.R = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.R
                                        g = Basic.Gibbs(aux.p, T)
                                        aux.rho = 1/g.gp
                                        aux.h = g.g-g.T*g.gT-Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.h_off
                                        aux.s = -g.gT-Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.s_off
                                        aux.cp = -g.T*g.gTT
                                        aux.cv = aux.cp-g.T/g.gpp*g.gTp^2
                                        aux.vt = 1/aux.rho*g.gTp/g.gp
                                        aux.vp = 1/aux.rho*g.gpp/g.gp
                                        aux.pd = -1/(aux.rho*aux.rho*aux.vp)
                                        aux.pt = -aux.vt/aux.vp
                    
                    return aux
                end

                "Density as function or pressure and temperature"
                function rho_props_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    # returns rho::Float64,          # Density
                    )
                    rho = aux.rho
                    
                    return rho
                end

                "Density as function or pressure and temperature"
                function rho_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns rho::Float64,          # Density
                    )
                    rho = rho_props_pT(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T))
                    
                    return rho
                end

                "Derivative function of rho_pT"
                function rho_pT_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    # returns rho_der::Float64,      # Derivative of density
                    )
                    rho_der = (-aux.rho*aux.rho*aux.vp)*p_der+(-aux.rho*aux.rho*aux.vt)*T_der
                    
                    return rho_der
                end

                "Specific enthalpy as function or pressure and temperature"
                function h_props_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    # returns h::Float64,            # Specific enthalpy
                    )
                    h = aux.h
                    
                    return h
                end

                "Specific enthalpy as function or pressure and temperature"
                function h_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns h::Float64,            # Specific enthalpy
                    )
                    h = h_props_pT(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T))
                    
                    return h
                end

                "Derivative function of h_pT"
                function h_pT_der(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    p_der::Float64,                  # Derivative of pressure
                    T_der::Float64,                  # Derivative of temperature
                    # returns h_der::Float64,        # Derivative of specific enthalpy
                    )
                    h_der = (1/aux.rho-aux.T*aux.vt)*p_der+aux.cp*T_der
                    
                    return h_der
                end

                "Specific entropy as function of pressure and temperature"
                function s_props_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    # returns s::Float64,            # Specific entropy
                    )
                    s = aux.s
                    
                    return s
                end

                "Temperature as function of pressure and temperature"
                function s_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns s::Float64,            # Specific entropy
                    )
                    s = s_props_pT(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T))
                    
                    return s
                end

                "Isothermal compressibility factor as function of pressure and temperature"
                function kappa_props_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    aux::Common.AuxiliaryProperties,   # Auxiliary record
                    # returns kappa::Float64,        # Isothermal compressibility factor
                    )
                    kappa = -aux.vp*aux.rho
                    
                    return kappa
                end

                "Isothermal compressibility factor as function of pressure and temperature"
                function kappa_pT(
                    #=
                    @extends Modelica.Icons.Function()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    # returns kappa::Float64,        # Isothermal compressibility factor
                    )
                    kappa = kappa_props_pT(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T))
                    
                    return kappa
                end
            end

            "Henry's law constant"
            function beta_H(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns beta_H::Float64,       # Henry's law constant
                )
                A::Array{Float64, 1},         
                B::Array{Float64, 1},         
                C::Array{Float64, 1},         
                psi::Array{Float64, 1},       
                beta::Array{Float64, 1},      
                Tr::Float64,                  
                tau::Float64,                 
                if 
                    ((T<273.15) || 
                    (T>Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p)))
                    beta_H = 0
                else
                    for k in 1:3
                beta[k] = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(T)*exp(A[k]/Tr+B[k]*tau^(0.355)/Tr+C[k]*Tr^(-0.41)*exp(tau))
                end
                
                    beta_H = 1/1.01325*
                    (psi[1]/beta[1]+psi[2]/beta[2]+psi[3]/beta[3])
                end
                
                return beta_H
            end

            "Enhancement factor as function of pressure and temperature"
            function f_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns f::Float64,            # Vapor-pressure enhancement factor
                )
#=

                
                function f_res(
                    #=
                    @extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction()
                    =#
                    p::Float64,                      # Pressure
                    T::Float64,                      # Temperature
                    )
                    x::Float64,                   
                    p_ws::Float64,                
                    kappa_T::Float64,             
                    v_ws::Float64,                
                    beta_H::Float64,              
                    psi_ws::Float64,              
                    baa::Float64,                 
                    baw::Float64,                 
                    bww::Float64,                 
                    caaa::Float64,                
                    caaw::Float64,                
                    caww::Float64,                
                    cwww::Float64,                
                    const R_bar::Float64 = Missing,
                    p_ws = if (T>=273.16); Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(T) else Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub(T) end
                                        if (p<p_ws)
                        kappa_T = 0
                    else
                        kappa_T = if (T>=273.16); Modelica.Media.Water.IF97_Utilities.kappa_pT(p, T) else Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.kappa_pT(p, T) end
                    end
                                        v_ws = if (T>=273.16); Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.molarMass/Modelica.Media.Water.IF97_Utilities.rho_pT(p, T) else Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Constants.MM/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT(p, T) end
                                        beta_H = Modelica.Media.Air.ReferenceMoistAir.Utilities.beta_H(p, T)
                                        baa = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Baa_dT(0, T)
                                        baw = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Baw_dT(0, T)
                                        bww = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Bww_dT(0, T)
                                        caaa = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Caaa_dT(0, T)
                                        caaw = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Caaw_dT(0, T)
                                        caww = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Caww_dT(0, T)
                                        cwww = Modelica.Media.Air.ReferenceMoistAir.Utilities.VirialCoefficients.Cwww_dT(0, T)
                                        y = 
                        ((1+kappa_T*p_ws)*(p-p_ws)-kappa_T*(p^2-p_ws^2)/2)/(R_bar*T)*v_ws+log(1-beta_H*(1-x*p_ws/p)*p)+(1-x*p_ws/p)^2*p/(R_bar*T)*baa-2*(1-x*p_ws/p)^2*p/(R_bar*T)*baw-(p-p_ws-(1-x*p_ws/p)^2*p)/(R_bar*T)*bww+(1-x*p_ws/p)^3*p^2/(R_bar*T)^2*caaa+3*(1-x*p_ws/p)^2*(1-2*(1-x*p_ws/p))*p^2/(2*(R_bar*T)^2)*caaw-3*(1-x*p_ws/p)^2*x*p_ws/p*p^2/(R_bar*T)^2*caww-((3-2*x*p_ws/p)*(x*p_ws/p)^2*p^2-p_ws^2)/(2*(R_bar*T)^2)*cwww-(1-x*p_ws/p)^2*(-2+3*x*p_ws/p)*x*p_ws/p*p^2/(R_bar*T)^2*baa*bww-2*(1-x*p_ws/p)^3*(-1+3*x*p_ws/p)*p^2/(R_bar*T)^2*baa*baw+6*(1-x*p_ws/p)^2*(x*p_ws/p)^2*p^2/(R_bar*T)^2*bww*baw-3*(1-x*p_ws/p)^4*p^2/(2*(R_bar*T)^2)*baa^2-2*(1-x*p_ws/p)^2*x*p_ws/p*(-2+3*x*p_ws/p)*p^2/(R_bar*T)^2*baw^2-(p_ws^2-(4-3*x*p_ws/p)*(x*p_ws/p)^3*p^2)/(2*(R_bar*T)^2)*bww^2-log(x)
                    
                    return 
                end
=#
                xmax::Float64,                
                if 
                    ((useEnhancementFactor==false) || 
                    (T>=Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p)))
                    f = 1
                else
                    xmax = if (T<273.16); 120.0 else 8.0 end
                    f = Modelica.Math.Nonlinear.solveOneNonlinearEquation()
                end
                
                return f
            end

            "Return density as a function of pressure p, temperature T and composition X"
            function rho_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns d::Float64,            # Density
                )
                pd::Float64,                  
                pl::Float64,                  
                xw::Float64,                  
                xws::Float64,                 
                if (X[1]==0)
                    d = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(p, T)
                else
                    xw = X[1]/(1-X[1])
                    xws = xws_pT(p, T)
                    pd = pd_pTX(p, T, X)
                    pl = p-pd
                    if ((xw<=xws) || (xws==-1))
                    if (T<273.16)
                    d = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+pd/
                    (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)
                else
                    d = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+IF97_new.rho_pT(pd, T)
                end
                else
                    if (T<273.16)
                    d = (1+xw)/
                    ((1+xws)/
                    (Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+pd/
                    (Modelica.Media.Air.ReferenceMoistAir.steam.R*T))+(xw-xws)/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT(p, T))
                else
                    d = (1+xw)/
                    ((1+xws)/
                    (Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+IF97_new.rho_pT(pd, T))+(xw-xws)/Modelica.Media.Water.IF97_Utilities.rho_pT(p, T))
                end
                end
                end
                
                return d
            end

            "Saturation partial pressure of steam"
            function pds_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns pds::Float64,          # Pressure
                )
                Tlim::Float64,                
                if (T>=273.16)
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat(T)
                    Tlim = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p)
                else
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub(T)
                    Tlim = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Tsub(p)
                end
                                if (T<=Tlim)
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.f_pT(p, T)*pds
                else
                    pds = -1
                end
                
                return pds
            end

            "partial pressure of steam"
            function pd_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns pd::Float64,           # partial pressure
                )
                xw::Float64,                  
                xws::Float64,                 
                pds::Float64,                 
                if (X[1]==0)
                    pd = 0
                else
                    xw = X[1]/(1-X[1])
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p, T)
                    pd = xw/
                    (Modelica.Media.Air.ReferenceMoistAir.k_mair+xw)*p
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    pd = if ((xw<=xws) || (xws==-1)); pd else pds end
                end
                
                return pd
            end

            "Humidity ration (absolute) of saturated humid air"
            function xws_pT(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                # returns xws::Float64,          # Absolute humidity ratio
                )
                pds::Float64,                 
                Tlim::Float64,                
                pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p, T)
                                Tlim = if (T<=273.16); Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Tsub(p) else Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p) end
                                xws = if (T<=Tlim); Modelica.Media.Air.ReferenceMoistAir.k_mair*pds/(p-pds) else -1 end
                
                return xws
            end

            "Relative humidity"
            function phi_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns phi::Float64,          # Relative humidity
                )
                xw::Float64,                  
                pd::Float64,                  
                pds::Float64,                 
                if (X[1]==0)
                    phi = 0
                else
                    xw = X[1]/(1-X[1])
                    if (T>=273.16)
                    pds = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T)
                else
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub(T)
                end
                    pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.f_pT(p, T)*pds
                    pd = xw/
                    (Modelica.Media.Air.ReferenceMoistAir.k_mair+xw)*p
                    if (pd<=pds)
                    phi = pd/pds
                else
                    phi = -1
                end
                end
                
                return phi
            end

            "Specific isobaric heat capacity"
            function cp_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns cp::Float64,           # Specific heat capacity
                )
                xw::Float64,                  
                xws::Float64,                 
                pd::Float64,                  
                pl::Float64,                  
                if (X[1]==0.0)
                    if (T>=773.15)
                    cp = Modelica.Media.Air.ReferenceAir.Air_Utilities.cp_pT(p, T)+Utilities.cp_dis_pTX(p, T, X)
                else
                    cp = Modelica.Media.Air.ReferenceAir.Air_Utilities.cp_pT(p, T)
                end
                else
                    pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                    pl = p-pd
                    xw = X[1]/(1-X[1])
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    if ((xw<=xws) || (xws==-1))
                    if (T>=773.15)
                    cp = X[1]*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cp_pT(pd, T)+X[2]*Modelica.Media.Air.ReferenceAir.Air_Utilities.cp_pT(pl, T)+Modelica.Media.Air.ReferenceMoistAir.Utilities.cp_dis_pTX(p, T, X)
                else
                    cp = X[1]*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cp_pT(pd, T)+X[2]*Modelica.Media.Air.ReferenceAir.Air_Utilities.cp_pT(pl, T)
                end
                else
                    cp = -1
                end
                end
                
                return cp
            end

            "Specific isochoric heat capacity"
            function cv_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns cv::Float64,           # Specific heat capacity
                )
                xw::Float64,                  
                xws::Float64,                 
                pd::Float64,                  
                pl::Float64,                  
                if (X[1]==0)
                    cv = Modelica.Media.Air.ReferenceAir.Air_Utilities.cv_pT(p, T)
                else
                    pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                    pl = p-pd
                    xw = X[1]/(1-X[1])
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    if ((xw<=xws) || (xws==-1))
                    cv = X[1]*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.cv_pT(pd, T)+X[2]*Modelica.Media.Air.ReferenceAir.Air_Utilities.cv_pT(pl, T)
                else
                    cv = -1
                end
                end
                
                return cv
            end

            "Specific enthalpy of moist air"
            function h_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns h::Float64,            # Specific enthalpy
                )
                xw::Float64,                  
                xws::Float64,                 
                pd::Float64,                  
                pl::Float64,                  
                if (X[1]==0)
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p, T)
                else
                    pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                    pl = p-pd
                    xw = X[1]/(1-X[1])
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    if ((xw<=xws) || (xws==-1))
                    if (T>=773.15)
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(1+xw)*Modelica.Media.Air.ReferenceMoistAir.Utilities.h_dis_pTX(p, T, X)
                else
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)
                end
                else
                    if (T<273.16)
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(xw-xws)*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT(p, T)
                else
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(xw-xws)*Modelica.Media.Water.IF97_Utilities.h_pT(p, T)
                end
                end
                    h = h/(1+xw)
                end
                
                return h
            end

            
            function h_dis_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns u::Float64,            # Reaction index
                )
                uges::Float64,                
                invMMX::Array{Float64, 1},       # Inverses of molar weights
                Mmix::Float64,                   # Molar mass of mixture
                massFraction::Array{MassFraction, 1},   # Mass fractions of components
                Y::Array{Float64, 1},            # Mole fractions of individual components (H2O, N2, O2, Ar) of moist air
                if (useDissociation==false)
                    u = 0
                else
                    massFraction = [X[1],X[2]*Xi_Air[1],X[2]*Xi_Air[2],X[2]*Xi_Air[3]]
                    for i in 1:4
                invMMX[i] = 1/MMX[i]
                end
                
                    Mmix = 1/(massFraction*invMMX)
                    for i in 1:4
                Y[i] = Mmix*massFraction[i]/MMX[i]
                end
                
                    uges = 1+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)
                    u = -T^2*
                    (Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V2(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[2]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V3(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[3]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V4(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[4]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V5(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[5]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V6(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[6])/uges*sum(massFraction[j]/MMX[j] for j in 1:4)
                end
                
                return u
            end

            "Specific entropy of moist air"
            function s_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns s::Float64,            # Specific entropy
                )
                xw::Float64,                  
                xws::Float64,                 
                pd::Float64,                  
                pl::Float64,                  
                if (X[1]==0)
                    s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(p, T)
                else
                    pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                    pl = p-pd
                    xw = X[1]/(1-X[1])
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    if ((xw<=xws) || (xws==-1))
                    if (T>=773.15)
                    s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.s_pT(pd, T)+(1+xw)*Modelica.Media.Air.ReferenceMoistAir.Utilities.s_dis_pTX(p, T, X)
                else
                    s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.s_pT(pd, T)
                end
                else
                    if (T<273.16)
                    s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.s_pT(pd, T)+(xw-xws)*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.s_pT(p, T)
                else
                    s = Modelica.Media.Air.ReferenceAir.Air_Utilities.s_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.s_pT(pd, T)+(xw-xws)*Modelica.Media.Water.IF97_Utilities.s_pT(p, T)
                end
                end
                    s = s/(1+xw)
                end
                
                return s
            end

            "Internal energy"
            function u_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns u::Float64,            # Specific entropy
                )
                u = Modelica.Media.Air.ReferenceMoistAir.Utilities.h_pTX(p, T, X)-p/Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX(p, T, X)
                
                return u
            end

            
            function cp_dis_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns u::Float64,            # Reaction index
                )
                uges::Float64,                
                invMMX::Array{Float64, 1},       # Inverses of molar weights
                Mmix::Float64,                   # Molar mass of mixture
                massFraction::Array{MassFraction, 1},   # Mass fractions of components
                Y::Array{Float64, 1},            # Mole fractions of individual components (H2O, N2, O2, Ar) of moist air
                if (useDissociation==false)
                    u = 0
                else
                    massFraction = [X[1],X[2]*Xi_Air[1],X[2]*Xi_Air[2],X[2]*Xi_Air[3]]
                    for i in 1:4
                invMMX[i] = 1/MMX[i]
                end
                
                    Mmix = 1/(massFraction*invMMX)
                    for i in 1:4
                Y[i] = Mmix*massFraction[i]/MMX[i]
                end
                
                    uges = 1+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)
                    u = 
                    (Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V2(T)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V3(T)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V4(T)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V5(T)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V6(T))/uges*sum(massFraction[j]/MMX[j] for j in 1:4)
                end
                
                return u
            end

            
            function s_dis_pTX(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                # returns u::Float64,            # Reaction index
                )
                uges::Float64,                
                invMMX::Array{Float64, 1},       # Inverses of molar weights
                Mmix::Float64,                   # Molar mass of mixture
                massFraction::Array{MassFraction, 1},   # Mass fractions of components
                Y::Array{Float64, 1},            # Mole fractions of individual components (H2O, N2, O2, Ar) of moist air
                if (useDissociation==false)
                    u = 0
                else
                    massFraction = [X[1],X[2]*Xi_Air[1],X[2]*Xi_Air[2],X[2]*Xi_Air[3]]
                    for i in 1:4
                invMMX[i] = 1/MMX[i]
                end
                
                    Mmix = 1/(massFraction*invMMX)
                    for i in 1:4
                Y[i] = Mmix*massFraction[i]/MMX[i]
                end
                
                    uges = 1+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)
                    u = -T*
                    (Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U2(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V2(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[2]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U3(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V3(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[3]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U4(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V4(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[4]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U5(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V5(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[5]+Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.U6(p, T, Y)*Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.V6(T)/Modelica.Media.Air.ReferenceMoistAir.Utilities.ReactionIndices.BB[6])/uges*sum(massFraction[j]/MMX[j] for j in 1:4)
                end
                
                return u
            end

            "Derivative of partial pressure of steam"
            function pd_pTX_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                X_der::Array{Float64, 1},        # Derivative of mass fractions
                # returns pd_der::Float64,       # Derivative of partial pressure
                )
                xw::Float64,                  
                xws::Float64,                 
                pds::Float64,                 
                if (X[1]==0)
                    pd_der = 0
                else
                    xw = X[1]/(1-X[1])
                    xw_der = (X_der[1]*(1-X[1])+X[1]*X_der[1])/(1-X[1])^2
                    pds_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT_der(p, T, p_der, T_der)
                    pd_der = 
                    (xw_der*
                    (Modelica.Media.Air.ReferenceMoistAir.k_mair+xw)-xw*xw_der)*p+xw/
                    (Modelica.Media.Air.ReferenceMoistAir.k_mair+xw)*p_der
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    pd_der = if ((xw<=xws) || (xws==-1)); pd_der else pds_der end
                end
                
                return pd_der
            end

            "Derivative of humidity ration (absolute) of saturated humid air"
            function xws_pT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                # returns xws_der::Float64,      # Derivative of absolute humidity ratio
                )
                pds::Float64,                 
                Tlim::Float64,                
                pds = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT(p, T)
                                pds_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.pds_pT_der(p, T, p_der, T_der)
                                Tlim = if (T<=273.16); Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Tsub(p) else Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p) end
                                if (T<=Tlim)
                    xws_der = Modelica.Media.Air.ReferenceMoistAir.k_mair*
                    ((pds_der*(p-pds)+pds*pds_der)/(p-pds)^2)
                else
                    xws_der = 0
                end
                
                return xws_der
            end

            "Derivative of saturation partial pressure of steam"
            function pds_pT_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                # returns pds_der::Float64,      # Derivative of pressure
                )
                Tlim::Float64,                
                if (T>=273.16)
                    pds_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.psat_der(T, T_der)
                    Tlim = Modelica.Media.Air.ReferenceMoistAir.Utilities.Water95_Utilities.Tsat(p)
                else
                    pds_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.psub_der(T, T_der)
                    Tlim = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.Basic.Tsub(p)
                end
                                if (T<=Tlim)
                    pds_der = pds_der
                else
                    pds_der = 0
                end
                
                return pds_der
            end

            "Derivative of density as a function of pressure p, temperature T and composition X"
            function rho_pTX_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                X_der::Array{Float64, 1},        # Derivative of mass fractions
                # returns d_der::Float64,        # Derivative of density
                )
                pd::Float64,                  
                pl::Float64,                  
                xw::Float64,                  
                xws::Float64,                 
                o::Array{Float64, 1},         
                if (X[1]==0)
                    d_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT_der(p, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(p, T), p_der, T_der)
                else
                    xw = X[1]/(1-X[1])
                    xw_der = (X_der[1])/(1-X[1])^2
                    xws = xws_pT(p, T)
                    xws_der = xws_pT_der(p, T, p_der, T_der)
                    pd = pd_pTX(p, T, X)
                    pd_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX_der(p, T, X, p_der, T_der, X_der)
                    pl = p-pd
                    pl_der = p_der-pd_der
                    if ((xw<=xws) || (xws==-1))
                    if (T<273.16)
                    d_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(pl, T), pl_der, T_der)+Modelica.Media.Air.ReferenceMoistAir.steam.R*(pd_der*T-pd*T_der)/
                    (Modelica.Media.Air.ReferenceMoistAir.steam.R*T)^2
                else
                    d_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(pl, T), pl_der, T_der)+IF97_new.rho_pT_der(pd, T, pd_der, T_der)
                end
                else
                    if (T<273.16)
                    o[1] = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)
                    o[2] = Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT(p, T)
                    o[3] = 
                    ((1+xws)/
                    (Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+pd/
                    (Modelica.Media.Air.ReferenceMoistAir.steam.R*T))+(xw-xws)/Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT(p, T))
                    o[4] = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(p, T), p_der, T_der)
                    o[5] = (xws_der*o[1]-(1+xws)*o[4])/o[1]^2+(pd_der*T-pd*T_der)/Modelica.Media.Air.ReferenceMoistAir.steam.R/T^2+
                    (xw_der*o[2]-xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT_der(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T), p_der, T_der))/o[2]^2-
                    (xws_der*o[2]-xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.rho_pT_der(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T), p_der, T_der))/o[2]^2
                    d_der = (xw_der*o[3]-(1+xw)*o[5])/o[3]^2
                else
                    o[1] = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+IF97_new.rho_pT(pd, T)
                    o[2] = Modelica.Media.Water.IF97_Utilities.rho_pT(p, T)
                    o[3] = 
                    ((1+xws)/
                    (Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT(pl, T)+IF97_new.rho_pT(pd, T))+(xw-xws)/Modelica.Media.Water.IF97_Utilities.rho_pT(p, T))
                    o[4] = Modelica.Media.Air.ReferenceAir.Air_Utilities.rho_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(p, T), p_der, T_der)+IF97_new.rho_pT_der(pd, T, pd_der, T_der)
                    o[5] = (xws_der*o[1]-(1+xws)*o[4])/o[1]^2+
                    (xw_der*o[2]-xw*Modelica.Media.Water.IF97_Utilities.rho_pT_der(p, T, Modelica.Media.Water.IF97_Utilities.waterBaseProp_pT(p, T), p_der, T_der))/o[2]^2-
                    (xws_der*o[2]-xws*Modelica.Media.Water.IF97_Utilities.rho_pT_der(p, T, Modelica.Media.Water.IF97_Utilities.waterBaseProp_pT(p, T), p_der, T_der))/o[2]^2
                    d_der = (xw_der*o[3]-(1+xw)*o[5])/o[3]^2
                end
                end
                end
                
                return d_der
            end

            
            function h_dis_pTX_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                X_der::Array{Float64, 1},        # Derivative of mass fractions
                # returns u_der::Float64,        # Derivative of reaction index
                )
                uges::Float64,                
                invMMX::Array{Float64, 1},       # Inverses of molar weights
                Mmix::Float64,                   # Molar mass of mixture
                Mmix_der::Float64,               # Derivative of molar mass of mixture
                massFraction::Array{MassFraction, 1},   # Mass fractions of components
                massFraction_der::Array{Float64, 1},   # Derivative of mass fractions of components
                Y::Array{Float64, 1},            # Mole fractions of individual components (H2O, N2, O2, Ar) of moist air
                Y_der::Array{Float64, 1},        # Derivative of mole fractions of individual components (H2O, N2, O2, Ar) of moist air
                if (useDissociation==false)
                    u_der = 0
                else
                    massFraction = [X[1],X[2]*Xi_Air[1],X[2]*Xi_Air[2],X[2]*Xi_Air[3]]
                    massFraction_der = [X_der[1],X_der[2]*Xi_Air[1],X_der[2]*Xi_Air[2],X_der[2]*Xi_Air[3]]
                    for i in 1:4
                invMMX[i] = 1/MMX[i]
                end
                
                    Mmix = 1/(massFraction*invMMX)
                    Mmix_der = -1/(massFraction*invMMX)^2*massFraction_der*invMMX
                    for i in 1:4
                Y[i] = Mmix*massFraction[i]/MMX[i]
                Y_der[i] = 
                    (Mmix_der*massFraction[i]+Mmix*massFraction_der[i])/MMX[i]
                end
                
                    uges = 1+Utilities.ReactionIndices.U2(p, T, Y)+Utilities.ReactionIndices.U3(p, T, Y)+Utilities.ReactionIndices.U4(p, T, Y)+Utilities.ReactionIndices.U5(p, T, Y)+Utilities.ReactionIndices.U6(p, T, Y)
                    uges_der = Utilities.ReactionIndices.U2_der(p, T, Y, p_der, T_der, Y_der)+Utilities.ReactionIndices.U3_der(p, T, Y, p_der, T_der, Y_der)+Utilities.ReactionIndices.U4_der(p, T, Y, p_der, T_der, Y_der)+Utilities.ReactionIndices.U5_der(p, T, Y, p_der, T_der, Y_der)+Utilities.ReactionIndices.U6_der(p, T, Y, p_der, T_der, Y_der)
                    o[1] = 
                    (Utilities.ReactionIndices.U2_der(p, T, Y, p_der, T_der, Y_der)*Utilities.ReactionIndices.V2(T)+Utilities.ReactionIndices.U2(p, T, Y)*Utilities.ReactionIndices.V2_der(T, T_der))/Utilities.ReactionIndices.BB[2]
                    o[2] = 
                    (Utilities.ReactionIndices.U3_der(p, T, Y, p_der, T_der, Y_der)*Utilities.ReactionIndices.V3(T)+Utilities.ReactionIndices.U3(p, T, Y)*Utilities.ReactionIndices.V3_der(T, T_der))/Utilities.ReactionIndices.BB[3]
                    o[3] = 
                    (Utilities.ReactionIndices.U4_der(p, T, Y, p_der, T_der, Y_der)*Utilities.ReactionIndices.V4(T)+Utilities.ReactionIndices.U4(p, T, Y)*Utilities.ReactionIndices.V4_der(T, T_der))/Utilities.ReactionIndices.BB[4]
                    o[4] = 
                    (Utilities.ReactionIndices.U5_der(p, T, Y, p_der, T_der, Y_der)*Utilities.ReactionIndices.V5(T)+Utilities.ReactionIndices.U5(p, T, Y)*Utilities.ReactionIndices.V5_der(T, T_der))/Utilities.ReactionIndices.BB[5]
                    o[5] = 
                    (Utilities.ReactionIndices.U6_der(p, T, Y, p_der, T_der, Y_der)*Utilities.ReactionIndices.V6(T)+Utilities.ReactionIndices.U6(p, T, Y)*Utilities.ReactionIndices.V6_der(T, T_der))/Utilities.ReactionIndices.BB[6]
                    l = 
                    (Utilities.ReactionIndices.U2(p, T, Y)*Utilities.ReactionIndices.V2(T)/Utilities.ReactionIndices.BB[2]+Utilities.ReactionIndices.U3(p, T, Y)*Utilities.ReactionIndices.V3(T)/Utilities.ReactionIndices.BB[3]+Utilities.ReactionIndices.U4(p, T, Y)*Utilities.ReactionIndices.V4(T)/Utilities.ReactionIndices.BB[4]+Utilities.ReactionIndices.U5(p, T, Y)*Utilities.ReactionIndices.V5(T)/Utilities.ReactionIndices.BB[5]+Utilities.ReactionIndices.U6(p, T, Y)*Utilities.ReactionIndices.V6(T)/Utilities.ReactionIndices.BB[6])
                    o[6] = -2*T*
                    (Utilities.ReactionIndices.U2(p, T, Y)*Utilities.ReactionIndices.V2(T)/Utilities.ReactionIndices.BB[2]+Utilities.ReactionIndices.U3(p, T, Y)*Utilities.ReactionIndices.V3(T)/Utilities.ReactionIndices.BB[3]+Utilities.ReactionIndices.U4(p, T, Y)*Utilities.ReactionIndices.V4(T)/Utilities.ReactionIndices.BB[4]+Utilities.ReactionIndices.U5(p, T, Y)*Utilities.ReactionIndices.V5(T)/Utilities.ReactionIndices.BB[5]+Utilities.ReactionIndices.U6(p, T, Y)*Utilities.ReactionIndices.V6(T)/Utilities.ReactionIndices.BB[6])-T^2*sum(o[1:5])
                    o[7] = uges_der*sum(massFraction[j]/MMX[j] for j in 1:4)+uges*sum(massFraction_der[j]/MMX[j] for j in 1:4)
                    u_der = 
                    (o[6]*
                    (uges*sum(massFraction[j]/MMX[j] for j in 1:4))-l*o[7])/
                    (uges*sum(massFraction[j]/MMX[j] for j in 1:4))^2
                end
                
                return u_der
            end

            "Derivative of specific enthalpy of moist air"
            function h_pTX_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                X_der::Array{Float64, 1},        # Derivative of mass fractions
                # returns h_der::Float64,        # Derivative of specific enthalpy
                )
                h::Float64,                   
                xw::Float64,                  
                xws::Float64,                 
                pd::Float64,                  
                pl::Float64,                  
                if (X[1]==0)
                    h_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT_der(p, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(p, T), p_der, T_der)
                else
                    pd = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX(p, T, X)
                    pd_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.pd_pTX_der(p, T, X, p_der, T_der, X_der)
                    pl = p-pd
                    pl_der = p_der-pd_der
                    xw = X[1]/(1-X[1])
                    xw_der = (X_der[1])/(1-X[1])^2
                    xws = Modelica.Media.Air.ReferenceMoistAir.Utilities.xws_pT(p, T)
                    xws_der = xws_pT_der(p, T, p_der, T_der)
                    if ((xw<=xws) || (xws==-1))
                    if (T>=773.15)
                    h_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(pl, T), pl_der, T_der)+xw_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT_der(pd, T, 0, pd_der, T_der)+xw_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.h_dis_pTX(p, T, X)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.h_dis_pTX_der(p, T, X, p_der, T_der, X_der)
                else
                    h_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(pl, T), pl_der, T_der)+xw_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT_der(pd, T, 0, pd_der, T_der)
                end
                else
                    if (T<273.16)
                    h_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(pl, T), pl_der, T_der)+xws_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT_der(pd, T, 0, pd_der, T_der)+xw_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT(p, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT_der(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T), p_der, T_der)-xws_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT(p, T)-xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT_der(p, T, Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.ice09BaseProp_pT(p, T), p_der, T_der)
                else
                    h_der = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT_der(pl, T, Modelica.Media.Air.ReferenceAir.Air_Utilities.airBaseProp_pT(p, T), pl_der, T_der)+xws_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT_der(pd, T, 0, pd_der, T_der)+xw_der*Modelica.Media.Water.IF97_Utilities.h_pT(p, T)+xw*Modelica.Media.Water.IF97_Utilities.h_pT_der(p, T, Modelica.Media.Water.IF97_Utilities.waterBaseProp_pT(p, T), p_der, T_der)-xws_der*Modelica.Media.Water.IF97_Utilities.h_pT(p, T)-xws*Modelica.Media.Water.IF97_Utilities.h_pT_der(p, T, Modelica.Media.Water.IF97_Utilities.waterBaseProp_pT(p, T), p_der, T_der)
                end
                end
                    if ((xw<=xws) || (xws==-1))
                    if (T>=773.15)
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(1+xw)*Modelica.Media.Air.ReferenceMoistAir.Utilities.h_dis_pTX(p, T, X)
                else
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xw*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)
                end
                else
                    if (T<273.16)
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(xw-xws)*Modelica.Media.Air.ReferenceMoistAir.Utilities.Ice09_Utilities.h_pT(p, T)
                else
                    h = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(pl, T)+xws*Modelica.Media.Air.ReferenceMoistAir.Utilities.IF97_new.h_pT(pd, T)+(xw-xws)*Modelica.Media.Water.IF97_Utilities.h_pT(p, T)
                end
                end
                    h_der = (h_der*(1+xw)-h*xw_der)/(1+xw)^2
                end
                
                return h_der
            end

            "Derivative of internal energy"
            function u_pTX_der(
                #=
                @extends Modelica.Icons.Function()
                =#
                p::Float64,                      # Pressure
                T::Float64,                      # Temperature
                X::Array{Float64, 1},            # Mass fractions
                p_der::Float64,                  # Derivative of pressure
                T_der::Float64,                  # Derivative of temperature
                X_der::Array{Float64, 1},        # Derivative of mass fractions
                # returns u_der::Float64,        # Derivative of specific entropy
                )
                u_der = Modelica.Media.Air.ReferenceMoistAir.Utilities.h_pTX_der(p, T, X, p_der, T_der, X_der)-
                    (p_der*Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX(p, T, X)-p*Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX_der(p, T, X, p_der, T_der, X_der))/Modelica.Media.Air.ReferenceMoistAir.Utilities.rho_pTX(p, T, X)^2
                
                return u_der
            end
        end
    end
end
