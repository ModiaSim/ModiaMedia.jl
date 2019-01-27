Directory FluidModels contains various variants to model Fluid components with Modia 
using ModiaMedia for the thermofluid-property modeling.

- Fluid1: Connector with (Medium,m_flow,r,p,h)
- Fluid2: Connector with (Medium,m_flow,r,state)
- Fluid3: Connector with (Medium,m_flow,r,state) and Medium in state

These models are for testing the variants and Modias symbolic transformation techniques.
Once satisfactory, it is planned to make a new Julia package (called e.g. ModiaFluid or ModiaStreams)
with the most promising variant.