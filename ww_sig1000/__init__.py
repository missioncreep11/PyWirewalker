"""Python port of the Wirewalker Nortek Signature ADCP velocity/turbulence toolbox.

Reference implementation: ../WW_Velocity_Processing_SWOT/*.m (velocity) and Devon
Northcott's ProcessSingleProfile.m (turbulence). Reads raw .ad2cp via
`from mhkit import dolfyn`. See the top-level README and ww_sig1000/validation/.
"""
