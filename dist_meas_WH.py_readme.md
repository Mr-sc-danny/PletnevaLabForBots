# dist_meas_WH.py

Source: groups\WH\dist_meas_WH.py

Purpose: Compares three distance relationships across the 2042, 2042f, 551, and 552 datasets.
Usage: Run with the configured MDTraj environment: python dist_meas_WH.py.
Inputs: The datasets list contains trajectory filenames, topology filenames, and atom-selection queries.
Calculation: For every dataset, computes two primary distances and an optional extra pair, then calculates mean, standard deviation, minimum, and maximum values.
Output: Writes newDist.txt, distance_vs_time_all.png, and one distance_vs_time_<dataset>.png per processed dataset.

The script is organized into explicit input, calculation, and output sections where applicable.
