# extract_spectra.py

Source: Agilent_UV_Extractor\extract_spectra.py

Purpose: Extracts raw UV-Vis spectra from Agilent/HP ChemStation binary `.SD` or `.KD` files.
Usage: `python extract_spectra.py sample.SD` or `python extract_spectra.py sample.KD`.
Inputs: One ChemStation `.SD` or `.KD` binary file containing wavelength and absorbance records.
Calculation: Parses wavelength axes, detects `CHPUVObject` scan records, reads absorbance arrays and sample metadata, de-duplicates column names, and writes tabular CSV data.
Output: A CSV beside the input file with `Wavelength_nm` and one column per spectrum or elapsed-time scan.

The Python scripts are organized into explicit input, calculation, and output sections where applicable.
Documentation follows the section-oriented style used by established Python repositories such as https://github.com/pandas-dev/pandas.
