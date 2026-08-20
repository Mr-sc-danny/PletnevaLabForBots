from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

ROOT = Path(r"E:\Desktop\research\MD\slit")
WORK_DIR = Path(__file__).resolve().parent

MODELS = {
    "1ayg": {
        "pdb_path": ROOT / "EPR_angle" / "1ayg.pdb",
        "heme_resname": "HEC",
        "heme_resseq": 81,
        "met_resseq": 59,
        "his_resseq": 14,
        "title": "1ayg",
    },
    "2042_60": {
        "pdb_path": ROOT / "EPR_angle" / "2042_60.pdb",
        "heme_resname": "HEM",
        "heme_resseq": 1,
        "met_resseq": 78,
        "his_resseq": 16,
        "title": "2042_60",
    },
    "c551_60": {
        "pdb_path": ROOT / "EPR_angle" / "c551_60.pdb",
        "heme_resname": "HEM",
        "heme_resseq": 83,
        "met_resseq": 61,
        "his_resseq": 16,
        "title": "c551_60",
    },
}


def parse_atoms(pdb_path: Path) -> dict[tuple[str, int, str], np.ndarray]:
    atoms: dict[tuple[str, int, str], np.ndarray] = {}
    with pdb_path.open() as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            resseq = int(line[22:26].strip())
            atoms[(resname, resseq, name)] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    return atoms


def fit_plane_normal(points: np.ndarray) -> np.ndarray:
    centroid = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - centroid)
    return vt[-1]


def project_onto_plane(vector: np.ndarray, normal: np.ndarray) -> np.ndarray:
    return vector - np.dot(vector, normal) * normal


def unit(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm < 1e-10:
        raise ValueError(f"Cannot normalize near-zero vector: {vector}")
    return vector / norm


def met_ligand_plane_projection(
    pdb_path: Path,
    heme_resname: str,
    heme_resseq: int,
    met_resseq: int,
    his_resseq: int,
) -> dict[str, np.ndarray | float]:
    atoms = parse_atoms(pdb_path)

    def get(resname: str, resseq: int, atom_name: str) -> np.ndarray:
        key = (resname, resseq, atom_name)
        if key not in atoms:
            raise KeyError(f"Atom not found in PDB: {key}")
        return atoms[key]

    fe = get(heme_resname, heme_resseq, "FE")
    na = get(heme_resname, heme_resseq, "NA")
    nb = get(heme_resname, heme_resseq, "NB")
    nc = get(heme_resname, heme_resseq, "NC")
    nd = get(heme_resname, heme_resseq, "ND")
    cg = get("MET", met_resseq, "CG")
    sd = get("MET", met_resseq, "SD")
    ce = get("MET", met_resseq, "CE")
    his_ne2 = get("HIS", his_resseq, "NE2")

    heme_normal = fit_plane_normal(np.array([na, nb, nc, nd]))
    if np.dot(heme_normal, sd - fe) < 0:
        heme_normal = -heme_normal

    proj_ch2 = unit(project_onto_plane(sd - cg, heme_normal))
    proj_ch3 = unit(project_onto_plane(sd - ce, heme_normal))
    bisector = unit(proj_ch2 + proj_ch3)
    ligand_trace = unit(np.cross(heme_normal, bisector))
    axis_na_nc = unit(project_onto_plane(nc - na, heme_normal))
    axis_nb_nd = unit(project_onto_plane(nd - nb, heme_normal))

    def plane_angle(vec: np.ndarray, ref: np.ndarray) -> float:
        cos_val = np.clip(abs(np.dot(vec, ref)), 0.0, 1.0)
        return float(np.degrees(np.arccos(cos_val)))

    return {
        "heme_normal": heme_normal,
        "proj_ch2": proj_ch2,
        "proj_ch3": proj_ch3,
        "bisector": bisector,
        "ligand_trace": ligand_trace,
        "angle_na_nc": plane_angle(ligand_trace, axis_na_nc),
        "angle_nb_nd": plane_angle(ligand_trace, axis_nb_nd),
        "fe_s_distance": float(np.linalg.norm(sd - fe)),
        "his_fe_distance": float(np.linalg.norm(his_ne2 - fe)),
    }


def write_report(model: str, result: dict[str, np.ndarray | float]) -> Path:
    report_path = WORK_DIR / f"epr_orientation_{model}.txt"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write(f"EPR Met axial-ligand orientation: {model}\n")
        handle.write("=" * 64 + "\n\n")
        handle.write(f"Heme plane normal: {result['heme_normal']}\n")
        handle.write(f"proj(CG→SD): {result['proj_ch2']}\n")
        handle.write(f"proj(CE→SD): {result['proj_ch3']}\n")
        handle.write(f"Ligand trace: {result['ligand_trace']}\n")
        handle.write(f"Fe–S distance: {result['fe_s_distance']:.3f} Å\n")
        handle.write(f"His–Fe distance: {result['his_fe_distance']:.3f} Å\n")
        handle.write(f"Angle vs NA→NC: {result['angle_na_nc']:.2f}°\n")
        handle.write(f"Angle vs NB→ND: {result['angle_nb_nd']:.2f}°\n")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze methionine ligand-plane orientation relative to a heme plane.")
    parser.add_argument("--model", choices=sorted(MODELS), default="1ayg")
    args = parser.parse_args()

    model = MODELS[args.model]
    result = met_ligand_plane_projection(
        pdb_path=model["pdb_path"],
        heme_resname=model["heme_resname"],
        heme_resseq=model["heme_resseq"],
        met_resseq=model["met_resseq"],
        his_resseq=model["his_resseq"],
    )
    report_path = write_report(args.model, result)
    print(f"Saved {report_path}")
    print(f"Angle vs NA→NC: {result['angle_na_nc']:.2f}°")
    print(f"Angle vs NB→ND: {result['angle_nb_nd']:.2f}°")


if __name__ == "__main__":
    main()
