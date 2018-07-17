## Welcome to QCL N2opocules

This is a tools for assessing N2O isotopocule fluxes and quality assurance from QCL data. It's developped by the Sustainable Agroecosystem (SAE) group at ETH Zurich. 

Scripts and sections:

- QCL-main.R
  - Define isotope standards (international and lab standards)
  - Read in the raw data from the QCL (implemented for Aerodyne QCLs - .stc and .str files) ... and merge
  - Aggregate the data to user defined timesteps
  - Define time vectors/formats
  - Calculate isotope ratios and delta values
  - Offset correction
  - Dilution correction
  - Calculate fluxes by chamber
  - Calculate source signature (Keeling plots)
  
- main_graphs.R:
  - Fig. 1: Dilution validation plot (checks how the dilution calibration was performed)
  - Fig. 2: Anchor-span plot (checks how the corrections perform against the reference gas)
  - Fig. 3: Keeling plot (Shows the concentration increase of a certain chamber with the corresponding keeling plot for the measured species)
  - Fig. 4.1: Chamber fluxes vs time
  - Fig. 4.2: Flux size vs. keeling plot delta values
- MonteMatti.R:
  - Simulates concentration increase and estimates the reliability of the source signuature calculation by Keeling plots (mixing model).
  - Estimates minimal required flux for reliable delta values by Monte Carlo Simulation

## Installation

its only a raw R script, that has not been turned into a function or package yet!


## Datasets

The sample dataset available is from a Aerodyne QCL used in a greenhouse plant variety experiment measuring soil gas emissions by automated static chambers. 


## References 

no accepted publications yet

visit www.sae.ethz.ch 
Github: https://github.com/pyroman1337/QCL-n2opocules 
Questions and feedback to: matti.barthel@usys.ethz.ch 
