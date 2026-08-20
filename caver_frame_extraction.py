from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

CONFIGS = {
    "2042": {
        "dcd_path": ROOT / "2042" / "2042_all.dcd",
        "top_path": ROOT / "2042" / "ionized.pdb",
        "output_dir": WORK_DIR / "caver_inputs_2042_with_heme_1000_3000",
    },
    "2042f_new": {
        "dcd_path": ROOT / "2042f_new" / "60nsp.dcd",
        "top_path": ROOT / "2042f_new" / "ionized.pdb",
        "output_dir": WORK_DIR / "caver_inputs_2042f_with_heme_1000_3000",
    },
}

START_FRAME = 1000
END_FRAME = 3000


def extract_frames(model: str) -> Path:
    config = CONFIGS[model]
    dcd_path = config["dcd_path"]
    top_path = config["top_path"]
    output_dir = config["output_dir"]

    if not dcd_path.exists():
        raise FileNotFoundError(f"Missing DCD file: {dcd_path}")
    if not top_path.exists():
        raise FileNotFoundError(f"Missing topology PDB: {top_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for pdb_file in output_dir.glob("*.pdb"):
        pdb_file.unlink()

    top = md.load(str(top_path)).topology
    selection = top.select("protein or resname HEMOX or resname HEM")
    if selection.size == 0:
        raise ValueError("Selection is empty: protein or heme residue selection")

    with md.open(str(dcd_path)) as dcd_file:
        num_frames = len(dcd_file)

    if END_FRAME > num_frames:
        raise ValueError(f"Requested end frame {END_FRAME}, but trajectory has only {num_frames} frames.")

    start_idx = START_FRAME - 1
    end_idx = END_FRAME - 1
    saved = 0
    frame_idx = 0

    for chunk in md.iterload(str(dcd_path), top=str(top_path), chunk=100):
        for local_i in range(chunk.n_frames):
            if frame_idx < start_idx:
                frame_idx += 1
                continue
            if frame_idx > end_idx:
                break

            saved += 1
            chunk[local_i].atom_slice(selection).save_pdb(str(output_dir / f"{saved}.pdb"))
            frame_idx += 1
        if frame_idx > end_idx:
            break

    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract CAVER input frames from the 2042 trajectories.")
    parser.add_argument("--model", choices=sorted(CONFIGS), default="2042")
    args = parser.parse_args()
    output_dir = extract_frames(args.model)
    print(f"Extracted frames to {output_dir}")


if __name__ == "__main__":
    main()
