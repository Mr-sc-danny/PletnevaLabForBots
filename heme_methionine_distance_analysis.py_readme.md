# heme_methionine_distance_analysis.py

Source: Slit_temp\heme_al\heme_methionine_distance_analysis.py

Purpose: Measures heme-to-ligand distances and angular statistics for the 2042, 2042f_new, and pac551 trajectories.
Usage: `python heme_methionine_distance_analysis.py --model pac551`
Inputs: Model-specific DCD/PDB paths and axial sulfur residue numbers in `CONFIGS`.
Calculation: Selects heme, sulfur, and histidine atoms; streams trajectory chunks; computes distances and angles; and filters frames from 20 to 60 ns.
Output: `distance_analysis_<model>.txt` with mean and standard-deviation summaries for each distance and angle.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
