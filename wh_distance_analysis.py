from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

CONFIGS = {
    "2042": {
        "dcd_path": ROOT / "2042" / "2042_all.dcd",
        "top_path": ROOT / "2042" / "ionized.pdb",
        "mode": "pair",
        "query_1": "resSeq 27 and name NE1",
        "query_2": "resSeq 16 and name ND1",
        "output_txt": "WH_2042.txt",
        "output_png": "distance_vs_time_2042.png",
        "title": "resSeq 27 NE1 to resSeq 16 ND1",
    },
    "2042f": {
        "dcd_path": ROOT / "2042f" / "2042fall.dcd",
        "top_path": ROOT / "2042f" / "ionized.pdb",
        "mode": "minimum_heavy",
        "query_1": "resSeq 16 and name ND1",
        "query_2": "resSeq 27 and not element H",
        "output_txt": "WH_2042f.txt",
        "output_png": "distance_vs_time_2042f.png",
        "title": "resSeq 16 ND1 to minimum heavy-atom distance in resSeq 27",
    },
}

DT_NS = 0.02


def run_analysis(model: str) -> Path:
    config = CONFIGS[model]
    top = md.load_topology(str(config["top_path"]))

    if config["mode"] == "pair":
        atom_1 = top.select(config["query_1"])
        atom_2 = top.select(config["query_2"])
        if len(atom_1) == 0 or len(atom_2) == 0:
            raise ValueError(f"Atom selection failed for {model}:\n{config['query_1']}\n{config['query_2']}")
        idx_1, idx_2 = atom_1[0], atom_2[0]
        pair_provider = lambda chunk: md.compute_distances(chunk, [[idx_1, idx_2]])[:, 0] * 10.0
    else:
        atom_1 = top.select(config["query_1"])
        heavy_atoms = top.select(config["query_2"])
        if len(atom_1) == 0 or len(heavy_atoms) == 0:
            raise ValueError(f"Atom selection failed for {model}:\n{config['query_1']}\n{config['query_2']}")
        idx_1 = atom_1[0]
        pairs = [[idx_1, heavy_atom] for heavy_atom in heavy_atoms]
        pair_provider = lambda chunk: np.min(md.compute_distances(chunk, pairs) * 10.0, axis=1)

    dist_all = []
    time_all = []
    frames_so_far = 0

    for chunk in md.iterload(str(config["dcd_path"]), top=str(config["top_path"]), chunk=1000):
        dist_all.append(pair_provider(chunk))
        time_all.append((np.arange(chunk.n_frames) + frames_so_far) * DT_NS)
        frames_so_far += chunk.n_frames

    distances = np.concatenate(dist_all)
    time = np.concatenate(time_all)
    mask = (time >= 20) & (time <= 60)
    distances_filtered = distances[mask]
    time_filtered = time[mask]
    if len(time_filtered) == 0:
        distances_filtered = distances
        time_filtered = time

    output_txt = WORK_DIR / config["output_txt"]
    with output_txt.open("w", encoding="utf-8") as handle:
        handle.write(f"Distance Analysis ({config['title']})\n")
        handle.write("=" * 50 + "\n\n")
        handle.write("Time window: 20-60 ns\n")
        handle.write(f"Number of frames: {len(distances_filtered)}\n")
        handle.write(f"Mean distance: {np.mean(distances_filtered):.3f} Å\n")
        handle.write(f"Std deviation: {np.std(distances_filtered):.3f} Å\n")
        handle.write(f"Min distance: {np.min(distances_filtered):.3f} Å\n")
        handle.write(f"Max distance: {np.max(distances_filtered):.3f} Å\n")

    plt.figure(figsize=(10, 4))
    plt.plot(time, distances, lw=1.5)
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance (Å)")
    plt.tight_layout()
    plt.savefig(WORK_DIR / config["output_png"], dpi=150)
    plt.close()

    return output_txt


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine WH distance-analysis variants into one script.")
    parser.add_argument("--model", choices=sorted(CONFIGS), default="2042")
    args = parser.parse_args()
    output_txt = run_analysis(args.model)
    print(f"Saved {output_txt}")


if __name__ == "__main__":
    main()
