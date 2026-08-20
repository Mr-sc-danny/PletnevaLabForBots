import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Edit these to your specific atoms (MDTraj selection syntax)
# Examples:
#   'resSeq 35 and name CA'      (PDB residue number 35 alpha carbon)
#   'resName HEM and name FE'    (heme iron; resName might be HEM/HEME)

dt_ns = 0.02  # 0.02 ns/frame

# Define datasets.
# Base queries compute (1->2) and (3->2).
# Optional extra pair computes an additional distance for selected datasets.
datasets = [
    ('2042', '2042.dcd', '2042.pdb', 'resSeq 27 and name NE1', 'resSeq 16 and name ND1', 'resSeq 25 and name N', ('resSeq 78 and name SD', 'resSeq 81 and name ND2')),
    ('2042f', '2042f.dcd', '2042f.pdb', 'resSeq 27 and name CE1', 'resSeq 16 and name ND1', 'resSeq 25 and name N', None),
    ('551', '551.dcd', '551.pdb', 'resSeq 27 and name CE1', 'resSeq 16 and name ND1', 'resSeq 25 and name N', ('resSeq 61 and name SD', 'resSeq 64 and name ND2')),
    ('552', 'ht552.dcd', 'ionized.pdb', 'resSeq 25 and name CE1', 'resSeq 14 and name ND1', 'resSeq 23 and name N', None),
]

results = {}
output_file = 'newDist.txt'

with open(output_file, 'w', encoding='utf-8') as f:
    f.write('Distance Analysis (20-60 ns)\n')
    f.write('=' * 60 + '\n\n')

for dataset_name, dcd_path, top_path, atom_query_1, atom_query_2, atom_query_3, extra_pair in datasets:
    print(f"\nProcessing {dataset_name}...")
    
    try:
        top = md.load_topology(top_path)
    except Exception as e:
        print(f"  Error loading topology for {dataset_name}: {e}")
        continue
    
    sel1 = top.select(atom_query_1)
    sel2 = top.select(atom_query_2)
    sel3 = top.select(atom_query_3)
    if len(sel1) == 0 or len(sel2) == 0 or len(sel3) == 0:
        print(f"  Atom selection failed for {dataset_name}. Check queries:\n    {atom_query_1}\n    {atom_query_2}\n    {atom_query_3}")
        continue
    
    idx1, idx2, idx3 = sel1[0], sel2[0], sel3[0]

    extra_idx_a = None
    extra_idx_b = None
    extra_label = None
    if extra_pair is not None:
        extra_query_a, extra_query_b = extra_pair
        extra_sel_a = top.select(extra_query_a)
        extra_sel_b = top.select(extra_query_b)
        if len(extra_sel_a) == 0 or len(extra_sel_b) == 0:
            print(
                f"  Extra atom selection failed for {dataset_name}. "
                f"Skipping extra pair:\n    {extra_query_a}\n    {extra_query_b}"
            )
        else:
            extra_idx_a, extra_idx_b = extra_sel_a[0], extra_sel_b[0]
            extra_label = (extra_query_a, extra_query_b)
    
    dist12_all = []
    dist32_all = []
    dist_extra_all = []
    time_all = []
    frames_so_far = 0
    
    try:
        for chunk in md.iterload(dcd_path, top=top_path, chunk=1000):
            d12 = md.compute_distances(chunk, [[idx1, idx2]])[:, 0] * 10.0  # nm -> Å
            d32 = md.compute_distances(chunk, [[idx3, idx2]])[:, 0] * 10.0  # nm -> Å
            dist12_all.append(d12)
            dist32_all.append(d32)
            if extra_idx_a is not None and extra_idx_b is not None:
                d_extra = md.compute_distances(chunk, [[extra_idx_a, extra_idx_b]])[:, 0] * 10.0
                dist_extra_all.append(d_extra)
            t = (np.arange(chunk.n_frames) + frames_so_far) * dt_ns
            time_all.append(t)
            frames_so_far += chunk.n_frames
    except Exception as e:
        print(f"  Error loading trajectory for {dataset_name}: {e}")
        continue

    if len(dist12_all) == 0 or len(dist32_all) == 0 or len(time_all) == 0:
        print(f"  No trajectory frames read for {dataset_name}.")
        continue
    
    # Concatenate
    d12 = np.concatenate(dist12_all)
    d32 = np.concatenate(dist32_all)
    d_extra = np.concatenate(dist_extra_all) if len(dist_extra_all) > 0 else None
    t = np.concatenate(time_all)
    
    # Filter data between 20-60 ns
    mask = (t >= 20) & (t <= 60)
    d12_filtered = d12[mask]
    d32_filtered = d32[mask]
    d_extra_filtered = d_extra[mask] if d_extra is not None else None
    t_filtered = t[mask]

    used_window = '20-60 ns'
    if len(t_filtered) == 0:
        d12_filtered = d12
        d32_filtered = d32
        if d_extra is not None:
            d_extra_filtered = d_extra
        t_filtered = t
        used_window = 'full trajectory (fallback: no frames in 20-60 ns)'

    def calc_stats(arr):
        return {
            'mean': np.mean(arr),
            'std': np.std(arr),
            'min': np.min(arr),
            'max': np.max(arr),
        }
    
    stats12 = calc_stats(d12_filtered)
    stats32 = calc_stats(d32_filtered)
    stats_extra = calc_stats(d_extra_filtered) if d_extra_filtered is not None else None
    
    results[dataset_name] = {
        'time': t,
        'distance': d32,
        't_filtered': t_filtered,
        'd12_filtered': d12_filtered,
        'd32_filtered': d32_filtered,
        'stats12': stats12,
        'stats32': stats32,
        'd_extra_filtered': d_extra_filtered,
        'stats_extra': stats_extra,
    }
    
    # Write all datasets into one file
    with open(output_file, 'a', encoding='utf-8') as f:
        f.write(f'Distance Analysis for {dataset_name}\n')
        f.write(f'Time window: {used_window}\n')
        f.write(f'Number of frames: {len(t_filtered)}\n\n')

        f.write(f'Atoms (1 -> 2): {atom_query_1} to {atom_query_2}\n')
        f.write(f'Mean distance: {stats12["mean"]:.3f} Å\n')
        f.write(f'Std deviation: {stats12["std"]:.3f} Å\n')
        f.write(f'Min distance: {stats12["min"]:.3f} Å\n')
        f.write(f'Max distance: {stats12["max"]:.3f} Å\n\n')

        f.write(f'Atoms (3 -> 2): {atom_query_3} to {atom_query_2}\n')
        f.write(f'Mean distance: {stats32["mean"]:.3f} Å\n')
        f.write(f'Std deviation: {stats32["std"]:.3f} Å\n')
        f.write(f'Min distance: {stats32["min"]:.3f} Å\n')
        f.write(f'Max distance: {stats32["max"]:.3f} Å\n')
        if stats_extra is not None and extra_label is not None:
            f.write('\n')
            f.write(f'Extra atoms: {extra_label[0]} to {extra_label[1]}\n')
            f.write(f'Mean distance: {stats_extra["mean"]:.3f} Å\n')
            f.write(f'Std deviation: {stats_extra["std"]:.3f} Å\n')
            f.write(f'Min distance: {stats_extra["min"]:.3f} Å\n')
            f.write(f'Max distance: {stats_extra["max"]:.3f} Å\n')
        f.write('-' * 60 + '\n\n')

    if stats_extra is not None:
        print(
            f"  {dataset_name} [{used_window}] - (1->2): {stats12['mean']:.3f} +- {stats12['std']:.3f} Å | "
            f"(3->2): {stats32['mean']:.3f} +- {stats32['std']:.3f} Å | "
            f"(extra): {stats_extra['mean']:.3f} +- {stats_extra['std']:.3f} Å"
        )
    else:
        print(
            f"  {dataset_name} [{used_window}] - (1->2): {stats12['mean']:.3f} +- {stats12['std']:.3f} Å | "
            f"(3->2): {stats32['mean']:.3f} +- {stats32['std']:.3f} Å"
        )

# Plot all datasets
if results:
    plt.figure(figsize=(12, 6))
    colors = {'2042': 'blue', '2042f': 'orange', '551': 'green'}
    for dataset_name, data in results.items():
        plt.plot(data['time'], data['distance'], lw=1, label=dataset_name, color=colors.get(dataset_name, 'gray'), alpha=0.7)
    
    plt.xlabel('Time (ns)')
    plt.ylabel('Distance (Å)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('distance_vs_time_all.png', dpi=150)
    print("\nSaved comparison plot: distance_vs_time_all.png")
    print(f"Saved distance summaries: {output_file}")
    plt.show()
    
    # Individual plots
    for dataset_name, data in results.items():
        plt.figure(figsize=(10, 4))
        plt.plot(data['time'], data['distance'], lw=1.5, color=colors.get(dataset_name, 'gray'))
        plt.xlabel('Time (ns)')
        plt.ylabel('Distance (Å)')
        plt.title(f'{dataset_name} - Distance vs Time')
        plt.tight_layout()
        plt.savefig(f'distance_vs_time_{dataset_name}.png', dpi=150)
        plt.close()
        print(f"Saved individual plot: distance_vs_time_{dataset_name}.png")