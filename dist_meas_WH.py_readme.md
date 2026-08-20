# dist_meas_WH.py

Source: Slit_temp\WH\dist_meas_WH.py

Purpose: Compares three related distance measurements across 2042, 2042f, 551, and 552 datasets.
Usage: Run `python dist_meas_WH.py` after setting the dataset paths and queries.
Inputs: The `datasets` list contains DCD paths, topology paths, and atom-selection expressions.
Calculation: Loads each topology, computes two primary distances and optional extra distances, filters 20-60 ns or falls back to the full trajectory, and calculates mean, standard deviation, minimum, and maximum.
Output: `newDist.txt`, `distance_vs_time_all.png`, and individual dataset distance plots.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
