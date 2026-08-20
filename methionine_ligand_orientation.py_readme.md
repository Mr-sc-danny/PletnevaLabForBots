# methionine_ligand_orientation.py

Source: groups\epr_angle\methionine_ligand_orientation.py

Purpose: Measures the orientation of an axial methionine ligand plane relative to the heme plane.
Usage: python methionine_ligand_orientation.py --model 2042_60
Inputs: A model PDB file and model-specific heme, methionine, and histidine residue settings in MODELS.
Calculation: Parses coordinates, fits the heme plane, projects C-S vectors, calculates the ligand-plane trace, and measures reference-axis angles.
Output: Writes epr_orientation_<model>.txt and prints the two orientation angles.

The script is organized into explicit input, calculation, and output sections where applicable.
