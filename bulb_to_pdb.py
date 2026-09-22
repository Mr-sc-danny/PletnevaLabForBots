"""
Insert CASTp bulb.json spheres as dummy atoms into an existing PDB.
Dummy atoms are inserted immediately before END.
Naming convention:
PBD submitted to Castp: protein.pdb
Bulb JSON: protein.json
Visualization comand:
size byattribute bfactor :DUM 
"""

import json

protein = "0019"

input_json = f"{protein}.json"
input_pdb = f"{protein}.pdb"
output_name = f"{protein}_with_bulb.pdb"


def build_dummy_atom_line(atom_num: int, sphere_index: int, sphere: dict) -> str:
    c = sphere["c"]
    r = sphere["r"]

    serial_str = f"{atom_num:5d}"
    atom_name = " X  "
    res_name = "DUM"
    chain_id = "D"
    res_seq = f"{(sphere_index // 100) + 1:4d}"

    x_str = f"{c['x']:8.3f}"
    y_str = f"{c['y']:8.3f}"
    z_str = f"{c['z']:8.3f}"
    occ_str = f"{1.00:6.2f}"
    bfac_str = f"{r:6.2f}"
    elem_str = " X"

    return (
        f"ATOM  {serial_str}  {atom_name}{res_name} {chain_id}{res_seq}    "
        f"{x_str}{y_str}{z_str}{occ_str}{bfac_str}          {elem_str}\n"
    )


with open(input_json, "r") as json_file:
    raw_data = json.load(json_file)

spheres = raw_data[0]

with open(input_pdb, "r") as pdb_file:
    pdb_lines = pdb_file.readlines()

insert_index = None
for idx, line in enumerate(pdb_lines):
    if line.strip() == "END":
        insert_index = idx
        break

if insert_index is None:
    raise ValueError("END not found in input PDB.")

max_serial = 0
for line in pdb_lines:
    if line.startswith(("ATOM  ", "HETATM")):
        serial_field = line[6:11].strip()
        if serial_field.isdigit():
            max_serial = max(max_serial, int(serial_field))

next_serial = max_serial + 1
remark_lines = [
    "REMARK   Dummy atoms from CASTp bulb.json pocket\n",
    "REMARK   Each atom represents a sphere; B-factor column = radius\n",
]
dummy_lines = []

for i, sphere in enumerate(spheres):
    dummy_lines.append(build_dummy_atom_line(next_serial, i, sphere))
    next_serial += 1

updated_lines = remark_lines + pdb_lines[:insert_index] + dummy_lines + pdb_lines[insert_index:]

with open(output_name, "w") as out:
    out.writelines(updated_lines)

    
