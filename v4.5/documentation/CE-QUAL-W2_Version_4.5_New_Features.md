# New Features: CE-QUAL-W2 Version 4.5

This version includes many new features and upgrades and is not file compatible with earlier versions because of many new variables in the control file. Differences between versions are shown in Part 5 of the User Manual. The new model includes the following features, developed by Portland State University, with support by ERDC Environmental Laboratory and Portland District of U.S. Army Corps of Engineers for the development of many of these enhancements):

1. Atmospheric deposition of any state variable - model user provides mass loading to waterbodies. There is now no need to specify a flow, temperature, and concentration file since only a time series of mass per area per time is used as an input.
2. Ability to specify directly output of flow balance file, N and P mass balance file, water level file
3. New generic constituent source: sediment release
4. New state variables: water age, N2, dissolved total gas pressure, bacteria, CH4, Fe(II), FeOOH, Mn(II), MnO2, H2S. Many of these before were only operative using sediment diagenesis as generic constituents or were recommended generic constituents.
5. New derived variables: turbidity (correlated to TSS), Secchi disk (based on light extinction), un-ionized ammonia (based on temperature, pH, and total ammonia), and Total Dissolved Gas (TDG).
6. Ammonia volatilization is computed based on unionized ammonia volatilization rate and is a derived variable. pH must also be active as a derived variable.
7. Implementation of variable algal settling velocity including buoyancy effects from Overman (2019). This allows for predicting the variable velocity of cyanobacteria allowing for rise and fall of the cells during the day and night.
8. Zooplankton settling is a new zooplankton parameter.
9. Implementation of algal toxin production based on Garstecki (2021)
10. Ability to generate lake contours easily (elevation vs time for temperature and dissolved oxygen)
11. Ability to generate river contours easily (model segment or distance along river vs time) for temperature and dissolved oxygen)
12. The C groups dissolved organic C (DOC) and particulate organic C (POC) with both labile and refractory groups were added as an alternative to organic matter groups. The option is available to use either C or organic matter. This was done earlier by Zhong in an earlier version of W2.
13. Sediment diagenesis updates. Many updates were made allowing for multiple vertical layers, simplified calculation, and much faster computation than before. Now when sediment diagenesis is turned ON, both zero order and first order sediment models are turned OFF internally. Also, the sediment diagenesis input file format has been updated and integrated into the Excel master sheet. This work was performed by Dr. Zhang at PSU.
14. Removed internal minimum reaeration coefficient value and allowed model user to set a minimum value if required. For example, for waterbodies designated as ‘LAKE’ with zero wind, the model will now predict zero reaeration unless a minimum value is set. The reaeration coefficient for the surface layer is now an output in the Time Series file.
15. Updates to auto-port selection
16. SYSTDG input files updated to csv input format.
17. Updates to particle settling and P adsorption onto inorganic SS. A bug was fixed where dissolved oxygen controlled P adsorption to particles rather than just P adsorption to Fe.
18. All external control files are included in an Excel master sheet for ease of writing out as csv files.

## Planned Enhancements

The following model enhancements are planned for CE-QUAL-W2 along with their status:

| # | Item | Description | Status |
|---|------|-------------|--------|
| 1 | Fish Fish Bioenergetics model | Tested and working in research code |
| 2 | Sediment transport | Simple model in research code |
| 3 | Smart particle tracking – Fish model | Tested and working in research code |
| 4 | Simultaneous water level solution | Currently, water surface is solved branch-by-branch. The new technique will involve solving all water surfaces for the system or waterbody simultaneously. | Tested and working in research code |
| 5 | W3 | 3D version of W2 | Tested and working in research code |
| 6 | Hyporheic flow algorithm | Groundwater-surface water interaction | Tested and working in research code. |
| 7 | Sediment channel bottom heating algorithm | Dynamic heat transfer between channel bottom and stream | Tested and working in research code. |
