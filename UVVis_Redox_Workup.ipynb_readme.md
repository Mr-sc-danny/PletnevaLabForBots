# UVVis_Redox_Workup.ipynb

Source: Echem_process\UVVis_Redox_Workup.ipynb

Purpose: Notebook pipeline for converting UV-Vis spectroelectrochemical spectra into fitted midpoint reduction potentials.
Usage: Open the notebook, set the `ALL PARAMETERS` cell, run once for the oxidative sweep and once for the reductive sweep.
Inputs: A folder of voltage-named `.SPC` spectra plus wavelength-picker settings, sweep label, and one/two-potential fit guesses.
Calculation: Parses spectra, performs baseline correction and wavelength selection, calculates reduced/oxidized fractions, fits a Boltzmann/Nernst model, and reports 95% confidence intervals.
Output: `<sweep_label>.xlsx` with equation, formatted parameters, numeric parameters, and fit data; intermediate `Full_data.xlsx` and `PotentialFraction.xlsx` files are also written.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
