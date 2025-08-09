# Frequently Asked Questions

- [How do I determine the appropriate water quality model to use?](#how-do-i-determine-the-appropriate-water-quality-model-to-use)
    - [CE-QUAL-W2: 2D Vertically Stratified, Longitudinally Varying, Laterally Averaged, Water Quality Model](#ce-qual-w2-2d-vertically-stratified-longitudinally-varying-laterally-averaged-water-quality-model)
    - [HEC-ResSim: 1D Vertically Stratified, Horizontally Averaged, Water Quality Model](#hec-ressim-1d-vertically-stratified-horizontally-averaged-water-quality-model)
    - [HEC-RAS: 1D Vertically Averaged (unstratified), Longitudinally Varying, Water Quality Model](#hec-ras-1d-vertically-averaged-unstratified-longitudinally-varying-water-quality-model)
    - [HEC-RAS: 2D Vertically Averaged (unstratified), Horizontally Varying, Water Quality Model](#hec-ras-2d-vertically-averaged-unstratified-horizontally-varying-water-quality-model)

## How do I determine the appropriate water quality model to use?

Determination of the best water quality model to use for a particular reservoir, river reach, etc., depends both on the water body being modeled and the objectives of the project.

### CE-QUAL-W2: 2D Vertically Stratified, Longitudinally Varying, Laterally Averaged, Water Quality Model
From the system perspective, if one or more of the following conditions applies to the waterbody being modeled, then CE-QUAL-W2 (2D hydrodynamic and water quality model) is the recommended model.

- The water body (lake, reservoir, or estuary) is stratified
- The water quality and other ecosystem processes in the water body are strongly affected by hydrodynamics and transport processes
- An accurate assessment of the processes in the water body are needed

Simpler models could be used depending on the system and the objectives of the project or management plan.

### HEC-ResSim: 1D Vertically Stratified, Horizontally Averaged, Water Quality Model

If the primary objective of the project is to meet downstream environmental objectives, and the details of in-reservoir water quality are not of concern, then the 1D vertically stratified capabilities in HEC-ResSim will accurately represent both stratified and unstratified conditions and meet the project objectives. HEC-ResSim will not provide the spatial information or accuracy that CE-QUAL-W2 is capable of, but it can be an excellent model for simulating water quality for identifying system-level ecosystem impacts of various water management alternatives.

### HEC-RAS: 1D Vertically Averaged (unstratified), Longitudinally Varying, Water Quality Model

For free-flowing river reaches in unmanaged watersheds or for unstratified fast-flowing river reaches between dams in regulated systems, then a 1D vertically averaged (unstratified) river hydraulics model like HEC-RAS is well suited to simulate the hydraulics of the flow from upstream to downstream and accurately computing the velocities and water depths along the river channel. Therefore, HEC-RAS (1D) is the recommended water quality modeling tool for these reaches. An added advantage of HEC-RAS is that RAS models exist for (almost) all of the rivers in the U.S. and many of the rivers around the world. These are often actively updated by H&H modeling teams. Therefore, HEC-RAS-1D leverages these existing models and the availability of expertise.

### HEC-RAS: 2D Vertically Averaged (unstratified), Horizontally Varying, Water Quality Model

For unstratified complex flows over the watershed, flood flows in the floodplain, braided river systems, weakly stratified estuaries, or other relatively shallow unstratified water bodies, two-dimensional (2D) HEC-RAS models represent the state of the art for simulating the flows, velocities, and water depths. Under the [Ecosystem Management and Restoration Research Program (EMRRP)](https://emrrp.el.erdc.dren.mil), water quality capabilities are being developed that leverage existing 2D HEC-RAS models. The Corps Library for Enviromental Analysis and Restoration of Watersheds (ClearWater), developed and maintained by ERDC-EL, is being extended under this project to develop a stand-alone water quality modeling tool ("ClearWater-Riverine") that will:

1. Read the outputs from computed 2D HEC-RAS models
2. Simulate 2D unstratified water quality processes (temperature and chemical kinetics) across the water surface
3. Output the results for visualization and analysis.
4. Modelers may then visualize and analyze the results using HEC-RAS and RAS Mapper as well as by a libary of Python programs.
