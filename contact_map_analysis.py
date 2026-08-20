from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

FOCUS_LABELS = [
    "",
    "35",
    "",
    "",
    "",
    "",
    "40",
    "",
    "",
    "",
    "45",
    "",
    "",
    "",
    "",
    "50",
    "",
    "",
    "",
    "",
    "55",
    "",
    "",
]

DATASETS = {
    "1353": {
        "dcd_path": ROOT / "1353" / "1353_all.dcd",
        "top_path": ROOT / "1353" / "ionized.psf",
        "frame_slice": slice(1499, 3000),
        "cutoff_nm": 1.0,
        "labels": [
            "1",
            "",
            "",
            "",
            "5",
            "",
            "",
            "",
            "",
            "10",
            "",
            "",
            "",
            "",
            "15",
            "",
            "",
            "",
            "",
            "20",
            "",
            "",
            "",
            "",
            "25",
            "",
            "",
            "",
            "",
            "30",
            "",
            "",
            "",
            "",
            "35",
            "",
            "",
            "",
            "",
            "40",
            "",
            "",
            "",
            "",
            "45",
            "",
            "",
            "",
            "",
            "50",
            "",
            "",
            "",
            "",
            "55",
            "",
            "",
            "",
            "",
            "60",
            "",
            "",
            "",
            "",
            "65",
            "",
            "",
            "",
            "",
            "70",
            "",
            "",
            "",
            "",
            "75",
            "",
            "",
            "",
            "",
            "80",
            "",
            "",
            "",
            "",
            "85",
            "",
            "",
            "",
            "",
            "90",
            "",
            "",
            "",
            "",
            "95",
            "",
            "",
            "",
            "",
            "100",
            "",
            "",
            "",
            "heme",
        ],
    },
    "2042": {
        "dcd_path": ROOT / "2042" / "2042_all.dcd",
        "top_path": ROOT / "2042" / "ionized.psf",
        "frame_slice": slice(999, 3000),
        "cutoff_nm": 0.5,
        "labels": [
            "-5",
            "",
            "",
            "",
            "",
            "1",
            "",
            "",
            "",
            "5",
            "",
            "",
            "",
            "",
            "10",
            "",
            "",
            "",
            "",
            "15",
            "",
            "",
            "",
            "",
            "20",
            "",
            "",
            "",
            "",
            "25",
            "",
            "",
            "",
            "",
            "30",
            "",
            "",
            "",
            "",
            "35",
            "",
            "",
            "",
            "",
            "40",
            "",
            "",
            "",
            "",
            "45",
            "",
            "",
            "",
            "",
            "50",
            "",
            "",
            "",
            "",
            "55",
            "",
            "",
            "",
            "",
            "60",
            "",
            "",
            "",
            "",
            "65",
            "",
            "",
            "",
            "",
            "70",
            "",
            "",
            "",
            "",
            "75",
            "",
            "",
            "",
            "",
            "80",
            "",
            "",
            "",
            "",
            "85",
            "",
            "",
            "",
            "",
            "90",
            "",
            "",
            "",
            "",
            "95",
            "",
            "",
            "",
            "",
            "100",
            "",
            "",
            "",
            "heme",
        ],
    },
}

FOCUS_SLICE = slice(33, 56)


def plot_heatmap(matrix: np.ndarray, labels: list[str], output_path: Path, title: str | None = None) -> None:
    plt.clf()
    fig = plt.figure(figsize=(28, 24))
    ax = fig.add_subplot(111)
    im = sns.heatmap(matrix, cmap="jet", ax=ax)
    ax.grid(which="both", alpha=0.5)
    plt.xticks(range(len(labels)), labels, rotation=45, size=22)
    plt.yticks(range(len(labels)), labels, rotation=45, size=22)
    cbar = im.collections[0].colorbar
    cbar.ax.tick_params(labelsize=22)
    ax.invert_yaxis()
    if title:
        ax.set_title(title)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_contact_map(dataset: str) -> Path:
    config = DATASETS[dataset]
    dcd_path = config["dcd_path"]
    top_path = config["top_path"]
    traj = md.load_dcd(str(dcd_path), str(top_path))
    traj = traj[config["frame_slice"]]
    top = traj.topology

    residue_number = np.arange(106)
    contact_maps = []
    contact_distances = []

    for i in range(len(residue_number)):
        if i % 5 == 0:
            print(f"Processing iteration {i}...")
        contact_map = []
        contact_distance = []
        for j in range(len(residue_number)):
            if i == j:
                contacts = 0
            else:
                dist = md.compute_contacts(traj, [[residue_number[i], residue_number[j]]])
                array = np.asarray(dist[0]).astype(float)
                distance = np.average(array)
                contact_distance.append(distance)
                contact = np.where(array < config["cutoff_nm"], 1, 0)
                contacts = np.average(contact)
            contact_map.append(contacts)
        contact_maps.append(contact_map)
        contact_distances.append(contact_distance)

    final_map = np.asarray(contact_maps).astype(float)
    final_distance = np.asarray(contact_distances).astype(float)

    map_path = WORK_DIR / f"contact_map_{dataset}.dat"
    dist_path = WORK_DIR / f"contact_distance_{dataset}.dat"
    np.savetxt(map_path, final_map, fmt="%1.3f")
    np.savetxt(dist_path, final_distance, fmt="%1.3f")
    plot_heatmap(final_map, config["labels"], WORK_DIR / f"contact_map_{dataset}.png")
    return map_path


def plot_focus_map(dataset: str, map_path: Path | None = None) -> None:
    if map_path is None:
        map_path = WORK_DIR / f"contact_map_{dataset}.dat"
    final_map = np.loadtxt(map_path)
    focus_map = final_map[FOCUS_SLICE, FOCUS_SLICE]
    plot_heatmap(focus_map, FOCUS_LABELS, WORK_DIR / f"contact_map_focus_{dataset}.png", title=f"{dataset} contact map focus")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute or plot curated contact maps.")
    parser.add_argument("--dataset", choices=sorted(DATASETS), default="2042")
    parser.add_argument("--mode", choices=["compute", "focus", "both"], default="both")
    parser.add_argument("--map-file", type=Path, default=None)
    args = parser.parse_args()

    if args.mode in {"compute", "both"}:
        map_path = compute_contact_map(args.dataset)
    else:
        map_path = args.map_file

    if args.mode in {"focus", "both"}:
        plot_focus_map(args.dataset, map_path)


if __name__ == "__main__":
    main()
