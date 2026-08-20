# methionine_ligand_orientation.py

Source: Slit_temp\epr_angle\methionine_ligand_orientation.py

Purpose: Determines the projected orientation of an axial methionine ligand plane relative to the heme plane.
Usage: `python methionine_ligand_orientation.py --model 2042_60`
Inputs: Model PDB files and residue mappings in `MODELS` for `1ayg`, `2042_60`, and `c551_60`.
Calculation: Parses ATOM/HETATM coordinates, fits the heme plane, projects the C-S vectors, computes the ligand-plane trace, and measures angles to the NA-NC and NB-ND axes.
Output: `epr_orientation_<model>.txt` containing vectors, Fe-S and His-Fe distances, and orientation angles.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
