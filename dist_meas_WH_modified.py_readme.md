# dist_meas_WH_modified.py

Source: Slit_temp\WH\dist_meas_WH_modified.py

Purpose: Compares three atom-distance measurements across 2042, 551, and 552 datasets.
Usage: Run `python dist_meas_WH_modified.py` after setting the dataset paths and queries.
Inputs: The `datasets` list contains four atom-selection expressions per dataset.
Calculation: Computes three distances per trajectory, applies the 20-60 ns window with a fallback, and calculates summary statistics.
Output: `newDist_modified.txt`, `distance_vs_time_modified.png`, and one three-panel plot per dataset.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
