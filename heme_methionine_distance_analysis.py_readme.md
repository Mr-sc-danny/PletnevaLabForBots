# heme_methionine_distance_analysis.py

Source: groups\heme_al\heme_methionine_distance_analysis.py

Purpose: Calculates heme-to-ligand distances and angle statistics for the 2042, 2042f_new, and pac551 trajectories.
Usage: python heme_methionine_distance_analysis.py --model pac551
Inputs: Model-specific DCD and topology paths plus the axial sulfur residue number in CONFIGS.
Calculation: Selects the heme, sulfur, and histidine atoms, iterates through the trajectory in chunks, calculates distances and angles, and filters the 20-60 ns window.
Output: Writes distance_analysis_<model>.txt with mean and standard-deviation summaries.

The script is organized into explicit input, calculation, and output sections where applicable.
