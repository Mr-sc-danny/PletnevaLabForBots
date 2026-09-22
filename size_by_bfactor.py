from chimerax.atomic import all_atoms

def size_by_bfactor(session, resname="DUM", scale=1.0):
    atoms = all_atoms(session)
    target = atoms.filter(atoms.residues.names == resname)
    for a in target:
        if a.bfactor is not None:
            a.radius = a.bfactor * scale
            a.draw_mode = a.SPHERE_STYLE  # switch out of stick mode so radius shows

size_by_bfactor(session)