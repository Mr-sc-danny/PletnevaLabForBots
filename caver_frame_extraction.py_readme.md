# caver_frame_extraction.py

Source: Slit_temp\frame_extraction\caver_frame_extraction.py

Purpose: Extracts trajectory frames 1000-3000 and keeps protein plus heme atoms for CAVER input.
Usage: `python caver_frame_extraction.py --model 2042` or `--model 2042f_new`.
Inputs: Model-specific DCD/PDB paths in `CONFIGS`; frame bounds are `START_FRAME` and `END_FRAME`.
Calculation: Validates files and frame count, selects protein/heme atoms, streams chunks, and saves each selected frame.
Output: Numbered PDB files in `caver_inputs_<model>_with_heme_1000_3000`.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
