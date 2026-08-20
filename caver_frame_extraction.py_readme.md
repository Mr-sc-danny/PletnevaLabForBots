# caver_frame_extraction.py

Source: groups\frame_extraction\caver_frame_extraction.py

Purpose: Extracts frames 1000-3000 from the 2042 or 2042f_new trajectories and keeps protein plus heme atoms for CAVER input.
Usage: python caver_frame_extraction.py --model 2042
Inputs: Model-specific DCD and PDB paths in CONFIGS; START_FRAME and END_FRAME define the extraction range.
Calculation: Validates files and frame count, loads the topology, selects protein/heme atoms, streams trajectory chunks, and writes each selected frame.
Output: Creates numbered PDB files in the configured caver_inputs_<model>_with_heme_1000_3000 directory.

The script is organized into explicit input, calculation, and output sections where applicable.
