# Cislunar Navigation System Analysis and Simulation

This repository houses analysis and simulation tools for Cislunar Navigation Systems (CNS). Tools are primarily developed in MATLAB and Python, which can be found in separate folders in the root directory.

## `./MATLAB/`

Tools developed in MATLAB are used primarily for analysis of candidate constellations -- frequently generated through GMAT or provided via SPICE kernels. Example tools include: evaluation of constellation metrics; generation of halo orbits; and optimization of satellite true anomaly, $f$, spacing within orbits.

## `./Python/`

This subdirectory houses Python tools developed primarily for navigation simulations, as well as the requisite GMAT scripts. These simulations intake data on satellite constellations generated through GMAT scripts and simulate CDS users navigating using measurements from the given constellation integrated into various filters (batch, Kalman, extended Kalman).

### GMAT

GMAT scripts are typically housed in `./python/gmat/`, while data products are dispersed to relevant directories where they are used.