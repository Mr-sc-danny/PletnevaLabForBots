from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

MODES = {
    "2way": {
        "output_png": "rmsd_sasa_2way.png",
        "title": "2042 vs Pa c551",
        "rmsd": [
            (ROOT / "551_All" / "2042_rmsd.dat", "#00FF00", "Slit_2042"),
            (ROOT / "551_All" / "pa551_rmsd.dat", "#FF00FF", "Pa c551"),
        ],
        "sasa": [
            (ROOT / "551_All" / "2042_SASA.dat", "#00FF00", "Slit_2042"),
            (ROOT / "551_All" / "pa551_sasa.dat", "#FF00FF", "Pa c551"),
        ],
    },
    "4way": {
        "output_png": "rmsd_sasa_4way.png",
        "title": "2042f, 2042, Pa c551, Ht c552",
        "rmsd": [
            (ROOT / "rmsd_sasa" / "2042f_rmsd.dat", "#335566", "2042f"),
            (ROOT / "rmsd_sasa" / "2042_rmsd.dat", "#779977", "2042"),
            (ROOT / "rmsd_sasa" / "pa551_rmsd.dat", "#aa4477", "Pa c551"),
            (ROOT / "rmsd_sasa" / "ht552_rmsd.dat", "#bb9955", "Ht c552"),
        ],
        "sasa": [
            (ROOT / "rmsd_sasa" / "2042f_sasaNew.dat", "#335566", "2042f"),
            (ROOT / "rmsd_sasa" / "2042_sasaNew.dat", "#779977", "2042"),
            (ROOT / "rmsd_sasa" / "pa551_sasaNew.dat", "#aa4477", "Pa c551"),
            (ROOT / "rmsd_sasa" / "ht552_sasaNew.dat", "#bb9955", "Ht c552"),
        ],
    },
}


def load_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, skiprows=1)
    return data[:, 0] / 50.0, data[:, 1]


def plot_mode(mode: str) -> Path:
    config = MODES[mode]
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Arial"
    plt.rcParams["mathtext.it"] = "Arial:italic"
    plt.rcParams["mathtext.bf"] = "Arial:bold"

    fig = plt.figure(figsize=(12, 8))

    ax_r = plt.subplot(211)
    for path, color, label in config["rmsd"]:
        time, values = load_curve(path)
        ax_r.plot(time, values, color=color, linewidth=1.5, label=label)
    ax_r.set_ylabel("RMSD (Å)", fontweight="bold", fontsize=18, fontname="Arial")
    ax_r.set_xlim([0, 60])
    ax_r.set_xticks(np.arange(0, 61, 20))
    ax_r.set_ylim([0, 3])
    ax_r.set_yticks([0, 1, 2])
    ax_r.tick_params(labelsize=18)
    ax_r.set_xticklabels([])
    for label in ax_r.get_xticklabels() + ax_r.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontname("Arial")
        label.set_fontsize(18)
    ax_r.legend(frameon=False)

    ax_s = plt.subplot(212)
    for path, color, label in config["sasa"]:
        time, values = load_curve(path)
        ax_s.plot(time, values, color=color, linewidth=1.5, label=label)
    ax_s.set_xlabel("Time (ns)", fontweight="bold", fontsize=18, fontname="Arial")
    ax_s.set_ylabel("Heme SASA (Å²)", fontweight="bold", fontsize=18, fontname="Arial")
    ax_s.set_xlim([0, 60])
    ax_s.set_xticks(np.arange(0, 61, 20))
    ax_s.set_ylim([0, 150])
    ax_s.set_yticks([0, 50, 100])
    ax_s.tick_params(labelsize=18)
    for label in ax_s.get_xticklabels() + ax_s.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontname("Arial")
        label.set_fontsize(18)
    ax_s.legend(frameon=False)
    ax_r.yaxis.set_label_coords(-0.06, 0.5)
    ax_s.yaxis.set_label_coords(-0.06, 0.5)
    plt.subplots_adjust(hspace=0)

    output_path = WORK_DIR / config["output_png"]
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine the RMSD/SASA comparison plots into one script.")
    parser.add_argument("--mode", choices=sorted(MODES), default="4way")
    args = parser.parse_args()
    output_path = plot_mode(args.mode)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
