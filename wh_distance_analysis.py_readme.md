# wh_distance_analysis.py

Source: groups\WH\wh_distance_analysis.py

Purpose: Measures the WH distance for either a direct atom pair or the minimum distance from one atom to heavy atoms in another residue.
Usage: python wh_distance_analysis.py --model 2042
Inputs: Model-specific trajectory paths and MDTraj atom-selection queries in CONFIGS.
Calculation: Loads the topology, selects atoms, computes distances in chunks, filters the 20-60 ns window with a full-trajectory fallback, and calculates summary statistics.
Output: Writes WH_<model>.txt and distance_vs_time_<model>.png.

The script is organized into explicit input, calculation, and output sections where applicable.
