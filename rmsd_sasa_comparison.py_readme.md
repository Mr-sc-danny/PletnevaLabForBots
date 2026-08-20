# rmsd_sasa_comparison.py

Source: Slit_temp\rmsd_sasa\rmsd_sasa_comparison.py

Purpose: Compares RMSD and heme SASA curves in two-way or four-way figures.
Usage: `python rmsd_sasa_comparison.py --mode 4way`.
Inputs: Data table paths, labels, and plot colors in `MODES`.
Calculation: Loads numeric tables, converts frame numbers to nanoseconds, and renders aligned RMSD and SASA panels.
Output: `rmsd_sasa_2way.png` or `rmsd_sasa_4way.png`.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
