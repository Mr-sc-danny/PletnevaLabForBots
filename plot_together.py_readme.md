# plot_together.py

Source: Echem_process\plot_together.py

Purpose: Plots oxidative and reductive electrochemical UV-Vis fitted curves together and summarizes midpoint potentials.
Usage: Put `oxi.xlsx` and `red.xlsx` beside the script, then run `python plot_together.py`.
Inputs: Excel files produced by `UVVis_Redox_Workup.ipynb`, with `parameters`, `parameters_numeric`, and `fit_data` sheets.
Calculation: Detects one- or two-potential models independently, reconstructs fitted fraction curves, plots raw points and fitted curves, and exports formatted midpoint values.
Output: `oxi_red_together.png` and `midpoint_potentials.xlsx`.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
