# License for this file: MIT
# Copyright 2018-2019, DLR Institute of System Dynamics and Control
#

"""
   module Fluid4_TestSimpleAir.jl

Date: Feb. 25, 2019

Test module to evaluate how fluid components can be modeled in Julia. The test models
in this module are for ModiaMedia SimpleAir. Since Modia does not yet support 
automatic state selection, the test models of different media are slightly different, 
because state=false has to be manually set at some components.


The components are modelled based on the new approach described in 

- Zimmer, Bender, Pollok (2018): Robust Modeling of Directed Thermofluid Flows in Complex Networks.
  2nd Japanese Modelica Conference, Tokyo, Proceedings, pp 39-48. Download:
  https://www.modelica.org/events/modelica2018japan/conference-proceedings/modelica-final-proceedings-2018-Japan.pdf


# Main Author of this module
```

      /|      Dirk Zimmer, Martin Otter
  ___/_|___   German Aerospace Center (DLR)
 /  /  /  /   Institute of System Dynamics and Control
/__/__/__/    https://www.dlr.de/sr/en/
   | /
   |/   
        
```  
"""
module Fluid4_TestSimpleAir

include("Components.jl")
using .Fluid4_Components

using  Modia
using  ModiaMedia
import ModiaMedia.ModiaMath  # included via ModiaMedia, to avoid requirement to add it in the standard environment


const MediumUsedInExamples = getMedium("SimpleAir")


@model TestLine begin
    fixedSource = FixedSource_pT(p0=2e5, T0=400.0, Medium=MediumUsedInExamples)
    fixedSink   = FixedSink(p0=1.5e5)
    pipe        = ShortPipe(l=1.0, A=0.01)

    @equations begin
        connect(fixedSource.outPlug, pipe.inPlug)
        connect(pipe.outPlug, fixedSink.inPlug)
    end
end

result = simulate(TestLine, 0.02; logSimulation=true)
ModiaMath.plot(result, ["pipe.m_flow", "pipe.outPlug.state.p", "pipe.outPlug.state.T", "pipe.outPlug.r", "pipe.eta"], heading="Fluid3_TestSimpleAir.jl: TestLine", figure=1)



@model TestJunction begin
    fixedSource = FixedSource_pT(p0=1.0e5, T0=400.0, Medium=MediumUsedInExamples)
    fixedSink   = FixedSink(  p0=0.8e5)
    pipe1       = ShortPipe(k=3e5, l=1.0, A=0.01)
    pipe2       = ShortPipe(       l=1.0, A=0.01)
    junction    = Junction()
    splitter    = Splitter()

    @equations begin
        connect(pipe2.outPlug      , junction.inPlugA)
        connect(pipe1.outPlug      , junction.inPlugB)
        connect(junction.outPlugC  , fixedSink.inPlug)
        connect(splitter.outPlugB  , pipe1.inPlug)
        connect(splitter.outPlugC  , pipe2.inPlug)
        connect(fixedSource.outPlug, splitter.inPlugA) 
    end
end

result = simulate(TestJunction, 0.02; logSimulation=true)
ModiaMath.plot(result, [("junction.inPlugA.m_flow", "junction.inPlugB.m_flow" , "junction.outPlugC.m_flow"),
                        ("splitter.inPlugA.m_flow", "splitter.outPlugB.m_flow", "splitter.outPlugC.m_flow"),
                        ("pipe1.outPlug.state.p", "pipe2.outPlug.state.p")], heading="Fluid3_TestSimpleAir.jl: TestJunction", figure=2)


@model TestLine2 begin
    fixedSource = FixedSource_pT(p0=2.0e5, T0=400.0, Medium=MediumUsedInExamples)
    fixedSink   = FixedSink(p0=1.0e5)
    pipe1       = ShortPipe(l=1.0, A=0.01)
    pipe2       = ShortPipe(l=1.0, A=0.01)
    volume      = ClosedVolume(p0=1.5e5, T0=500.0, Medium=MediumUsedInExamples)

    @equations begin
        connect(fixedSource.outPlug, pipe1.inPlug)
        connect(pipe1.outPlug      , volume.inPlug)
        connect(volume.outPlug     , pipe2.inPlug)
        connect(pipe2.outPlug      , fixedSink.inPlug)
    end
end

result = simulate(TestLine2, 1.0, logTranslation=true, logSimulation=true, tearing=true)
ModiaMath.plot(result, [("volume.inPlug.m_flow", "volume.outPlug.m_flow"),
                         "volume.medium.p",
                        ("fixedSource.outPlug.state.T", "volume.medium.T"),
                        ("volume.h", "volume.medium.h"), 
                        ("volume.outPlug.state.d"),
                        ("volume.outPlug.state.eta")], 
                         heading="Fluid3_TestSimpleAir.jl: TestLine2", figure=3)

#=
@model TestLine3 begin
    fixedSource = FixedSource_pT(p0=2.0e5, T0=400.0, Medium=MediumUsedInExamples)
    fixedSink   = FixedSink(p0=1.0e5)
    pipe1       = ShortPipe(l=1.0, A=0.01)
    pipe2       = ShortPipe(l=1.0, A=0.01)
    volume1     = ClosedVolume(p0=1.5e5, T0=500.0, V=0.5, Medium=MediumUsedInExamples)
    volume2     = ClosedVolume(V=0.5, Medium=MediumUsedInExamples, 
                               medium = BaseProperties(MediumUsedInExamples, p_start=StaticPressure(start=1e5), T_start=Temperature(start=500.0)))

    @equations begin
        connect(fixedSource.outPlug, pipe1.inPlug)
        connect(pipe1.outPlug      , volume1.inPlug)
        connect(volume1.outPlug    , volume2.inPlug)        
        connect(volume2.outPlug    , pipe2.inPlug)
        connect(pipe2.outPlug      , fixedSink.inPlug)
    end
end

result = simulate(TestLine3, 1.0, logTranslation=true, logSimulation=true, tearing=true)
add!(result, "pipe2.outPlug", ["T", "d"])
add!(result, "fixedSource.outPlug", ["T", "d"])
ModiaMath.plot(result, [("volume1.inPlug.m_flow", "volume1.outPlug.m_flow"),
                         "volume1.medium.p",
                        ("fixedSource.outPlug.T", "volume1.medium.T"),
                        ("volume1.h", "volume1.medium.h")], 
                         heading="Fluid3_TestSimpleAir.jl: TestLine3", figure=4)
=#



@model TestCompressor begin
    fixedSource = FixedSource_pT(p0=1.0e5, T0=400, Medium=MediumUsedInExamples)
    fixedSink   = FixedSink(  p0=1.5e5)
    pipe1       = ShortPipe(l=1.0, A=0.01)
    pipe2       = ShortPipe(l=1.0, A=0.01)
    compressor  = Compressor(l=1.0, A=0.01)

    @equations begin
        connect(fixedSource.outPlug, pipe1.inPlug)
        connect(pipe1.outPlug      , compressor.inPlug)
        connect(compressor.outPlug , pipe2.inPlug)
        connect(pipe2.outPlug      , fixedSink.inPlug)
    end
end
println("\n\nSimulate TestCompressor")
# result = simulate(TestCompressor, 0.02; logTranslation=true, logSimulation=true, tearing=false)
result = simulate(TestCompressor, 0.02; logTranslation=true, logSimulation=true, tearing=true, aliasElimination=true, removeSingularities=false, automaticStateSelection=true)
ModiaMath.plot(result, [ "compressor.inPlug.m_flow",
                        ("compressor.inPlug.state.p", "compressor.outPlug.state.p", "pipe2.outPlug.state.p")], 
                         heading="Fluid3_TestSimpleAir.jl: TestCompressor", figure=5)


@model TestLoop begin
    pipe1       = ShortPipe(       l=1.0, A=0.01)
    pipe2       = ShortPipe(k=3e5, l=1.0, A=0.01)
    compressor  = Compressor(l=1.0, A=0.01)
    volume      = ClosedVolume(p0=1e5, T0=500, Medium=MediumUsedInExamples)

    @equations begin
        connect(volume.outPlug    , pipe1.inPlug)
        connect(pipe1.outPlug     , compressor.inPlug)
        connect(compressor.outPlug, pipe2.inPlug)
        connect(pipe2.outPlug     , volume.inPlug)
    end
end

result = simulate(TestLoop, 0.02; logTranslation=true, logSimulation=true, tearing=false, aliasElimination=true, removeSingularities=false, automaticStateSelection=true)
ModiaMath.plot(result, [ "compressor.inPlug.m_flow",
                        ("compressor.inPlug.state.p", "compressor.outPlug.state.p", "pipe2.outPlug.state.p", "volume.medium.p")], 
                         heading="Fluid3_TestSimpleAir.jl: TestLoop", figure=6)


end
 