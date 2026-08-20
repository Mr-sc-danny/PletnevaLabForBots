import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Edit these to your specific atoms (MDTraj selection syntax)
# Examples:
#   'resSeq 35 and name CA'      (PDB residue number 35 alpha carbon)
#   'resName HEM and name FE'    (heme iron; resName might be HEM/HEME)

dt_ns = 0.02  # 0.02 ns/frame

# Define datasets.
# Computes distances from atom_query_1 to atom_query_2, atom_query_3, and atom_query_4.
datasets = [
    ('2042', '2042.dcd', '2042.pdb', 'resSeq 27 and name NE1', 'resSeq 16 and name NE2', 'resSeq 16 and name ND1', 'resSeq 12 and name O'),
    ('551', '551.dcd', '551.pdb', 'resSeq 27 and name OH', 'resSeq 16 and name NE2', 'resSeq 16 and name ND1', 'resSeq 12 and name O'),
    ('552', 'ht552.dcd', 'ionized.pdb', 'resSeq 25 and name OH', 'resSeq 14 and name NE2', 'resSeq 14 and name ND1', 'resSeq 10 and name O'),
]

results = {}
output_file = 'newDist_modified.txt'

with open(output_file, 'w', encoding='utf-8') as f:
    f.write('Distance Analysis (20-60 ns)\n')
    f.write('=' * 60 + '\n\n')

for dataset_name, dcd_path, top_path, atom_query_1, atom_query_2, atom_query_3, atom_query_4 in datasets:
    print(f"\nProcessing {dataset_name}...")
    
    try:
        top = md.load_topology(top_path)
    except Exception as e:
        print(f"  Error loading topology for {dataset_name}: {e}")
        continue
    
    sel1 = top.select(atom_query_1)
    sel2 = top.select(atom_query_2)
    sel3 = top.select(atom_query_3)
    sel4 = top.select(atom_query_4)
    if len(sel1) == 0 or len(sel2) == 0 or len(sel3) == 0 or len(sel4) == 0:
        print(f"  Atom selection failed for {dataset_name}. Check queries:\n    {atom_query_1}\n    {atom_query_2}\n    {atom_query_3}\n    {atom_query_4}")
        continue
    
    idx1, idx2, idx3, idx4 = sel1[0], sel2[0], sel3[0], sel4[0]
    
    dist12_all = []
    dist13_all = []
    dist14_all = []
    time_all = []
    frames_so_far = 0
    
    try:
        for chunk in md.iterload(dcd_path, top=top_path, chunk=1000):
            d12 = md.compute_distances(chunk, [[idx1, idx2]])[:, 0] * 10.0  # nm -> Å
            d13 = md.compute_distances(chunk, [[idx1, idx3]])[:, 0] * 10.0  # nm -> Å
            d14 = md.compute_distances(chunk, [[idx1, idx4]])[:, 0] * 10.0  # nm -> Å
            dist12_all.append(d12)
            dist13_all.append(d13)
            dist14_all.append(d14)
            t = (np.arange(chunk.n_frames) + frames_so_far) * dt_ns
            time_all.append(t)
            frames_so_far += chunk.n_frames
    except Exception as e:
        print(f"  Error loading trajectory for {dataset_name}: {e}")
        continue

    if len(dist12_all) == 0 or len(dist13_all) == 0 or len(dist14_all) == 0 or len(time_all) == 0:
        print(f"  No trajectory frames read for {dataset_name}.")
        continue
    
    # Concatenate
    d12 = np.concatenate(dist12_all)
    d13 = np.concatenate(dist13_all)
    d14 = np.concatenate(dist14_all)
    t = np.concatenate(time_all)
    
    # Filter data between 20-60 ns
    mask = (t >= 20) & (t <= 60)
    d12_filtered = d12[mask]
    d13_filtered = d13[mask]
    d14_filtered = d14[mask]
    t_filtered = t[mask]

    used_window = '20-60 ns'
    if len(t_filtered) == 0:
        d12_filtered = d12
        d13_filtered = d13
        d14_filtered = d14
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
    stats13 = calc_stats(d13_filtered)
    stats14 = calc_stats(d14_filtered)
    
    results[dataset_name] = {
        'time': t,
        't_filtered': t_filtered,
        'd12_filtered': d12_filtered,
        'd13_filtered': d13_filtered,
        'd14_filtered': d14_filtered,
        'stats12': stats12,
        'stats13': stats13,
        'stats14': stats14,
    }
    
    # Write all datasets into one file
    with open(output_file, 'a', encoding='utf-8') as f:
        f.write(f'Distance Analysis for {dataset_name}\n')
        f.write(f'Time window: {used_window}\n')
        f.write(f'Number of frames: {len(t_filtered)}\n\n')

        f.write(f'Distance 1: {atom_query_1} to {atom_query_2}\n')
        f.write(f'Mean distance: {stats12["mean"]:.3f} Å\n')
        f.write(f'Std deviation: {stats12["std"]:.3f} Å\n')
        f.write(f'Min distance: {stats12["min"]:.3f} Å\n')
        f.write(f'Max distance: {stats12["max"]:.3f} Å\n\n')

        f.write(f'Distance 2: {atom_query_1} to {atom_query_3}\n')
        f.write(f'Mean distance: {stats13["mean"]:.3f} Å\n')
        f.write(f'Std deviation: {stats13["std"]:.3f} Å\n')
        f.write(f'Min distance: {stats13["min"]:.3f} Å\n')
        f.write(f'Max distance: {stats13["max"]:.3f} Å\n\n')

        f.write(f'Distance 3: {atom_query_1} to {atom_query_4}\n')
        f.write(f'Mean distance: {stats14["mean"]:.3f} Å\n')
        f.write(f'Std deviation: {stats14["std"]:.3f} Å\n')
        f.write(f'Min distance: {stats14["min"]:.3f} Å\n')
        f.write(f'Max distance: {stats14["max"]:.3f} Å\n')
        f.write('-' * 60 + '\n\n')

    print(
        f"  {dataset_name} [{used_window}] - Distance 1: {stats12['mean']:.3f} +- {stats12['std']:.3f} Å | "
        f"Distance 2: {stats13['mean']:.3f} +- {stats13['std']:.3f} Å | "
        f"Distance 3: {stats14['mean']:.3f} +- {stats14['std']:.3f} Å"
    )

# Plot all datasets
if results:
    colors = {'2042': 'blue', '551': 'green', '552': 'red'}
    labels = {
        '2042': 'NE2 (blue), ND1 (cyan), O (lightblue)',
        '551': 'NE2 (green), ND1 (lime), O (lightgreen)',
        '552': 'NE2 (red), ND1 (pink), O (lightcoral)',
    }
    
    # Combined plot for all three distances
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    for dataset_name, data in results.items():
        color = colors.get(dataset_name, 'gray')
        axes[0].plot(data['t_filtered'], data['d12_filtered'], lw=1, label=dataset_name, color=color, alpha=0.7)
        axes[1].plot(data['t_filtered'], data['d13_filtered'], lw=1, label=dataset_name, color=color, alpha=0.7)
        axes[2].plot(data['t_filtered'], data['d14_filtered'], lw=1, label=dataset_name, color=color, alpha=0.7)
    
    axes[0].set_ylabel('Distance 1 (Å)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].set_ylabel('Distance 2 (Å)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    axes[2].set_xlabel('Time (ns)')
    axes[2].set_ylabel('Distance 3 (Å)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('distance_vs_time_modified.png', dpi=150)
    print("\nSaved comparison plot: distance_vs_time_modified.png")
    print(f"Saved distance summaries: {output_file}")
    plt.show()
    
    # Individual plots per dataset
    for dataset_name, data in results.items():
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))
        color = colors.get(dataset_name, 'gray')
        
        axes[0].plot(data['t_filtered'], data['d12_filtered'], lw=1.5, color=color)
        axes[0].set_ylabel('Distance 1 (Å)')
        axes[0].set_title(f'{dataset_name} - Distance 1 vs Time')
        axes[0].grid(True, alpha=0.3)
        
        axes[1].plot(data['t_filtered'], data['d13_filtered'], lw=1.5, color=color)
        axes[1].set_ylabel('Distance 2 (Å)')
        axes[1].set_title(f'{dataset_name} - Distance 2 vs Time')
        axes[1].grid(True, alpha=0.3)
        
        axes[2].plot(data['t_filtered'], data['d14_filtered'], lw=1.5, color=color)
        axes[2].set_xlabel('Time (ns)')
        axes[2].set_ylabel('Distance 3 (Å)')
        axes[2].set_title(f'{dataset_name} - Distance 3 vs Time')
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'distance_vs_time_modified_{dataset_name}.png', dpi=150)
        plt.close()
        print(f"Saved individual plot: distance_vs_time_modified_{dataset_name}.png")
