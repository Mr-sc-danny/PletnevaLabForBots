from pathlib import Path
import mdtraj as md

dcd_path = Path(r\"E:\Desktop\research\MD\slit\mouth\ht552\ht552_all.dcd\")
top_path = Path(r\"E:\Desktop\research\MD\slit\mouth\ht552\ionized.pdb\")
out_dir = Path(r\"E:\Desktop\research\MD\slit\mouth\ht552\caver_inputs_ht552_from_dcd\")
start_frame = 1
end_frame = 2
stride = 1

if not dcd_path.exists():
    raise FileNotFoundError(f\"Missing DCD file: {dcd_path}\")
if not top_path.exists():
    raise FileNotFoundError(f\"Missing topology PDB: {top_path}\")

for old_pdb in out_dir.glob("*.pdb"):
    old_pdb.unlink()

top = md.load(str(top_path)).topology
selection = top.select("protein or resname HEMOX or resname HEM")
if selection.size == 0:
    raise ValueError("Selection is empty: protein or heme residue selection")

with md.open(str(dcd_path)) as dcd_file:
    total_frames = len(dcd_file)

if end_frame == 0:
    end_frame = total_frames

if end_frame > total_frames:
    raise ValueError(f\"Requested end frame {end_frame}, but trajectory has only {total_frames} frames.\")

start_idx = start_frame - 1
end_idx = end_frame - 1

saved = 0
frame_idx = 0
for chunk in md.iterload(str(dcd_path), top=str(top_path), chunk=100):
    for local_i in range(chunk.n_frames):
        if frame_idx < start_idx:
            frame_idx += 1
            continue
        if frame_idx > end_idx:
            break
        if ((frame_idx - start_idx) % stride) != 0:
            frame_idx += 1
            continue

        saved += 1
        out_path = out_dir / f\"{saved}.pdb\"
        chunk[local_i].atom_slice(selection).save_pdb(str(out_path))
        frame_idx += 1

    if frame_idx > end_idx:
        break

print(f\"Extracted {saved} frames from {dcd_path.name} to {out_dir}\")
