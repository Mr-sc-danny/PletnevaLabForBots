"""
Plot Oxidative and Reductive fitted curves together, with their raw data points.

Refactored from Echem_plot_together.m.

Usage: place this script in the same folder as `oxi.xlsx` and `red.xlsx` (produced by
UVVis_Redox_Workup.ipynb, one run with sweep_label = "oxi" and one with sweep_label = "red"),
then run:

    python plot_together.py

Edit the parameters block below if your filenames or display range differ.

Each of oxi.xlsx / red.xlsx is read independently and can be either a one-potential fit
(A1, Ema, n1) or a two-potential fit (adds A2, Emb, n2) -- whichever was used when that file
was produced. The script detects this per file and prints which model it found.

Also writes `midpoint_potentials.xlsx`: one row per sweep (oxi/red), with Ema (and Emb if
either file is a two-potential fit) reported as "value \u00b1 95% CI error", matching the
notebook's own reporting exactly.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =================================================================================
# PARAMETERS -- edit as needed
# =================================================================================
oxi_filename = "oxi.xlsx"      # output of the notebook with sweep_label = "oxi"
red_filename = "red.xlsx"      # output of the notebook with sweep_label = "red"

voltage_range = (-400, 700)    # range the fitted curves are drawn over (mV)
VoltageDisplay = (-250, 600)   # x-axis display window (mV)
# =================================================================================

CONST = 39.585 / 1000   # n*F/(R*T) at room temperature, in 1/mV -- same constant as in the notebook


def load_fit_file(path):
    """Read one sweep's fitted parameters + the Potential/Fraction data used to fit it.

    Expects the file layout written by the notebook's Part 3: a 'parameters_numeric'
    sheet (parameter, value) and a 'fit_data' sheet (Potential, Fraction).
    """
    params_df = pd.read_excel(path, sheet_name="parameters_numeric")
    params = dict(zip(params_df["parameter"], params_df["value"]))
    data_df = pd.read_excel(path, sheet_name="fit_data")
    return params, data_df["Potential"].to_numpy(), data_df["Fraction"].to_numpy()


def get_param_display(path, name):
    """Read the formatted 'value \u00b1 95% CI error' string for one parameter from the
    notebook's 'parameters' sheet (same sheet/format as its own fit_results reporting,
    so the number shown here always matches what the notebook printed).

    Returns None if that parameter isn't in this file (e.g. Emb for a one-potential fit).
    """
    params_df = pd.read_excel(path, sheet_name="parameters")
    match = params_df[params_df["parameter"] == name]
    if match.empty:
        return None
    return match.iloc[0]["value \u00b1 95% CI error"]


def describe_model(params):
    """Human-readable label for whether this file is a one- or two-potential fit,
    purely for the console message printed in main() -- doesn't affect the plot."""
    if "A2" in params:
        return "two-potential (A1, Ema, n1, A2, Emb, n2)"
    return "one-potential (A1, Ema, n1)"


def make_fraction_func(params):
    """Build Fraction(V) from whatever parameters are present. Each sweep is checked
    independently, so oxi.xlsx and red.xlsx can be a mix of one- and two-potential fits."""
    def fraction(V):
        total = params["A1"] / (np.exp((V - params["Ema"]) * params["n1"] * CONST) + 1)
        if "A2" in params:
            total = total + params["A2"] / (np.exp((V - params["Emb"]) * params["n2"] * CONST) + 1)
        return total
    return fraction


def save_potential_summary(path="midpoint_potentials.xlsx"):
    """Write a small side-by-side table of the midpoint potential(s) with their 95% CI
    error, exactly as reported in oxi.xlsx/red.xlsx -- one row per sweep, with an Emb
    column added only if at least one sweep used a two-potential fit (left blank for
    whichever sweep only has one potential)."""
    oxi_ema = get_param_display(oxi_filename, "Ema")
    red_ema = get_param_display(red_filename, "Ema")
    oxi_emb = get_param_display(oxi_filename, "Emb")
    red_emb = get_param_display(red_filename, "Emb")

    rows = [
        {"": "oxi", "Ema": oxi_ema},
        {"": "red", "Ema": red_ema},
    ]
    if oxi_emb is not None or red_emb is not None:
        rows[0]["Emb"] = oxi_emb if oxi_emb is not None else ""
        rows[1]["Emb"] = red_emb if red_emb is not None else ""

    summary_df = pd.DataFrame(rows)

    # Overwrite the file even if it already exists; fall back to a differently-named
    # file only if it's currently open elsewhere (e.g. in Excel) and locked.
    try:
        summary_df.to_excel(path, index=False)
        print(f"Saved {path}")
    except PermissionError:
        fallback = path.replace(".xlsx", "_new.xlsx")
        summary_df.to_excel(fallback, index=False)
        print(f"'{path}' is open elsewhere (e.g. in Excel) and could not be overwritten -- "
              f"saved to '{fallback}' instead.")


def main():
    # Load both sweeps and report which equation (one- or two-potential) each used.
    oxi_params, oxi_V, oxi_F = load_fit_file(oxi_filename)
    red_params, red_V, red_F = load_fit_file(red_filename)

    print(f"{oxi_filename}: detected {describe_model(oxi_params)}")
    print(f"{red_filename}: detected {describe_model(red_params)}")

    save_potential_summary()

    # Reconstruct each sweep's fitted curve as a callable function of voltage.
    oxi_fraction = make_fraction_func(oxi_params)
    red_fraction = make_fraction_func(red_params)

    v_plot = np.linspace(voltage_range[0], voltage_range[1], 400)

    fig, ax = plt.subplots(figsize=(7, 5))

    # Fitted curves (fplot equivalent) + raw data points for both sweeps, matching
    # the original MATLAB script's combined plot.
    ax.plot(v_plot, oxi_fraction(v_plot), linewidth=1.5, label="Oxidative")
    ax.plot(v_plot, red_fraction(v_plot), linewidth=1.5, label="Reductive")
    ax.scatter(oxi_V, oxi_F, s=15, color="black")
    ax.scatter(red_V, red_F, s=15, color="black")

    ax.set_xlim(VoltageDisplay)
    ax.set_xlabel("Voltage (mV)", fontweight="bold")
    ax.set_ylabel("Fraction of Fe$^{II}$", fontweight="bold")
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight("bold")

    ax.legend(loc="upper right", prop={"weight": "bold"})

    plt.tight_layout()
    plt.savefig("oxi_red_together.png", dpi=200)   # saved since a plain script has no inline display
    print("Saved oxi_red_together.png")
    plt.show()


if __name__ == "__main__":
    main()