# rmsd_sasa_comparison.py

Source: groups\rmsd_sasa\rmsd_sasa_comparison.py

Purpose: Creates paired RMSD and heme-SASA comparison plots for two-way or four-way dataset comparisons.
Usage: python rmsd_sasa_comparison.py --mode 4way
Inputs: The MODES configuration identifies the RMSD/SASA data tables, colors, and labels.
Calculation: Loads data tables, converts frames to nanoseconds, plots RMSD and SASA panels, and applies consistent axis styling.
Output: Writes rmsd_sasa_2way.png or rmsd_sasa_4way.png.

The script is organized into explicit input, calculation, and output sections where applicable.
