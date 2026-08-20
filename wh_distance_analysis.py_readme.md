# wh_distance_analysis.py

Source: Slit_temp\WH\wh_distance_analysis.py

Purpose: Runs the WH distance analysis in either direct pair mode or minimum-heavy-atom mode.
Usage: `python wh_distance_analysis.py --model 2042` or `python wh_distance_analysis.py --model 2042f`.
Inputs: Trajectory paths and MDTraj atom-selection queries in `CONFIGS`.
Calculation: Selects atoms, calculates distances in chunks, filters the 20-60 ns interval with a full-trajectory fallback, and summarizes the distribution.
Output: `WH_<model>.txt` and `distance_vs_time_<model>.png`.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
