# contact_map_analysis.py

Source: Slit_temp\contact_maps\contact_map_analysis.py

Purpose: Calculates residue-residue contact frequencies and mean pair distances from an MD trajectory, then renders full and focused heatmaps.
Usage: `python contact_map_analysis.py --dataset 2042 --mode both`
Inputs: DCD and PSF paths plus dataset-specific frame ranges, cutoffs, and labels in `DATASETS`.
Calculation: Loads the trajectory, computes every residue-pair distance, applies the configured nanometer cutoff, averages contact indicators, and optionally extracts the focused residue window.
Output: `contact_map_<dataset>.dat`, `contact_distance_<dataset>.dat`, `contact_map_<dataset>.png`, and `contact_map_focus_<dataset>.png`.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
