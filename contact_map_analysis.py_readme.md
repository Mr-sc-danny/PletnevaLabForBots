# contact_map_analysis.py

Source: groups\contact_maps\contact_map_analysis.py

Purpose: Computes residue-residue contact frequencies and average distances from MD trajectories, then creates full and focused heatmaps.
Usage: python contact_map_analysis.py --dataset 2042 --mode both
Inputs: Dataset paths and residue-label configuration are defined in DATASETS; supported datasets are 1353 and 2042.
Calculation: Loads the trajectory, computes pairwise contact frequencies using the configured cutoff, and optionally extracts the focused residue window.
Output: Writes contact_map_<dataset>.dat, contact_distance_<dataset>.dat, contact_map_<dataset>.png, and contact_map_focus_<dataset>.png.

The script is organized into explicit input, calculation, and output sections where applicable.
