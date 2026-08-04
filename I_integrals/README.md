# Electron Radial I-Integral Calculation

This directory contains the numerical code used to calculate the electron
radial integrals \(I^{(1)}\) and \(I^{(2)}\) for Schwarzschild black-hole
wave functions.

The calculation reads precomputed electron wave functions, reflection and
transmission coefficients, evaluates the radial integrals for different
electron energies and angular quantum numbers, and stores the results in
FITS files.

## Files

- `electron01.sh`  
  SLURM array-job submission script. It divides the \(h'\) grid into
  computation blocks and starts a separate array task for each block.
  Array task 0 waits for the calculation tasks to finish and then merges
  the output files.

- `makeintegrand.py`  
  Main calculation script. For a specified pair of grid indices `x` and
  `y`, it loads the electron wave functions and calculates the
  \(I^{(1)}\) and \(I^{(2)}\) integrals for
  \(\kappa=-10,\ldots,-1,1,\ldots,10\).

- `I_functions_class.py`  
  Contains the electron and photon radial wave-function classes, RK4
  integrators, reflection and transmission coefficient calculations, and
  the numerical definitions of \(I^{(1)}\) and \(I^{(2)}\).

- `module1.py`  
  Provides the conversion between the Schwarzschild radial coordinate
  \(r\) and the tortoise coordinate \(r_*\).

- `merge_fits.py`  
  Sorts the individual output files by the \(h'\) index and combines them
  into one multi-extension FITS file.

- `Constants1.fits`  
  Stores physical and numerical parameters, including the black-hole mass,
  electron mass, Hawking temperature, integration tolerance, and paths to
  the precomputed electron and photon wave-function files.

- `plot.ipynb`  
  Notebook used to load the calculated integrals, construct correction
  terms, and produce diagnostic plots.

## Requirements

The code requires Python 3 and the following packages:

```bash
numpy
scipy
matplotlib
sympy
astropy
fitsio
