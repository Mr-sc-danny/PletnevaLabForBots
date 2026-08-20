# dist_meas_WH_modified.py

Source: groups\WH\dist_meas_WH_modified.py

Purpose: Compares three atom-distance measurements for the 2042, 551, and 552 datasets.
Usage: Run with the configured MDTraj environment: python dist_meas_WH_modified.py.
Inputs: The datasets list contains trajectory filenames, topology filenames, and four atom-selection queries.
Calculation: Computes three distances per trajectory, applies the 20-60 ns window with a full-trajectory fallback, and calculates summary statistics.
Output: Writes newDist_modified.txt, distance_vs_time_modified.png, and distance_vs_time_modified_<dataset>.png files.

The script is organized into explicit input, calculation, and output sections where applicable.
