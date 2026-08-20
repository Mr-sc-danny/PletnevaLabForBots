from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import numpy as np

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

CONFIGS = {
    "2042": {
        "dcd_path": ROOT / "2042" / "2042_all.dcd",
        "top_path": ROOT / "2042" / "ionized.pdb",
        "sd_resseq": 78,
        "output_name": "distance_analysis_2042.txt",
    },
    "2042f_new": {
        "dcd_path": ROOT / "2042f_new" / "60nsp.dcd",
        "top_path": ROOT / "2042f_new" / "ionized.pdb",
        "sd_resseq": 78,
        "output_name": "distance_analysis_2042f_new.txt",
    },
    "pac551": {
        "dcd_path": ROOT / "pac551" / "Pa_c551_ox_1.dcd",
        "top_path": ROOT / "pac551" / "ionized.pdb",
        "sd_resseq": 61,
        "output_name": "distance_analysis_pac551.txt",
    },
}

DT_NS = 0.02


def run_analysis(model: str) -> Path:
    config = CONFIGS[model]
    top = md.load_topology(str(config["top_path"]))
    queries = [
        "resname HEMO and name FE",
        "resname HEMO and name NA",
        "resname HEMO and name NB",
        "resname HEMO and name NC",
        "resname HEMO and name ND",
        f"resSeq {config['sd_resseq']} and name SD",
        "resSeq 16 and name NE2",
    ]
    selections = [top.select(query) for query in queries]
    if any(len(selection) == 0 for selection in selections):
        missing = [query for query, selection in zip(queries, selections) if len(selection) == 0]
        raise ValueError("Atom selection failed. Check queries:\n" + "\n".join(missing))

    idx = [selection[0] for selection in selections]
    pairs = [[idx[0], idx[i]] for i in range(1, 7)]
    angle_triplets = [[idx[6], idx[0], idx[i]] for i in range(1, 6)]

    dist_all = []
    angle_all = []
    time_all = []
    frames_so_far = 0

    for chunk in md.iterload(str(config["dcd_path"]), top=str(config["top_path"]), chunk=1000):
        dist_all.append(md.compute_distances(chunk, pairs) * 10.0)
        angle_all.append(md.compute_angles(chunk, angle_triplets) * (180.0 / np.pi))
        time_all.append((np.arange(chunk.n_frames) + frames_so_far) * DT_NS)
        frames_so_far += chunk.n_frames

    distances = np.concatenate(dist_all, axis=0)
    angles = np.concatenate(angle_all, axis=0)
    time = np.concatenate(time_all)
    mask = (time >= 20) & (time <= 60)
    distances_filtered = distances[mask]
    angles_filtered = angles[mask]

    output_path = WORK_DIR / config["output_name"]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("Distance Analysis (Atom 1 to Atoms 2-6)\n")
        handle.write("=" * 50 + "\n\n")
        handle.write("Time window: 20-60 ns\n")
        handle.write(f"Number of frames: {distances_filtered.shape[0]}\n\n")

        labels = [
            "FE to NA",
            "FE to NB",
            "FE to NC",
            "FE to ND",
            f"FE to SD (resid {config['sd_resseq']})",
            "FE to NE2 (resid 16)",
        ]
        for index, label in enumerate(labels):
            handle.write(f"{label}\n")
            handle.write(f"  Mean distance: {np.mean(distances_filtered[:, index]):.3f} Å\n")
            handle.write(f"  Std deviation: {np.std(distances_filtered[:, index]):.3f} Å\n\n")

        handle.write("Angle Analysis (7–1–2..6)\n")
        handle.write("=" * 50 + "\n\n")
        for index, label in enumerate(["Angle 7-1-2", "Angle 7-1-3", "Angle 7-1-4", "Angle 7-1-5", "Angle 7-1-6"]):
            handle.write(f"{label}\n")
            handle.write(f"  Mean angle: {np.mean(angles_filtered[:, index]):.3f} deg\n")
            handle.write(f"  Std deviation: {np.std(angles_filtered[:, index]):.3f} deg\n\n")

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine the heme-ligand distance analyses into one script.")
    parser.add_argument("--model", choices=sorted(CONFIGS), default="2042")
    args = parser.parse_args()
    output_path = run_analysis(args.model)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
