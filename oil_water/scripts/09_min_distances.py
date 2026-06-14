"""
Minimum pairwise distance analysis.

Computes per-frame minimum distances between:
  - CHI ↔ NAP  (closest atom pair between chitosan and any naphthalene)
  - CHI ↔ PAR  (closest atom pair between chitosan and any paraffin)

Input:  complex.prmtop, traj_centered_50ns.dcd
Output: min_distances.csv  (frame, time_ps, dmin_chi_nap, dmin_chi_par)
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

# ===============================
# Load
# ===============================
u = mda.Universe("complex.prmtop", "traj_centered_50ns.dcd")

sel_chi = u.select_atoms("resname CHI")
sel_nap = u.select_atoms("resname NAP")
sel_par = u.select_atoms("resname PAR")

# ===============================
# Per-frame minimum distances
# ===============================
records = []

for ts in u.trajectory:
    # CHI–NAP minimum distance
    pairs_nap, dists_nap = capped_distance(
        sel_chi.positions, sel_nap.positions,
        max_cutoff=100.0, box=ts.dimensions
    )
    dmin_nap = dists_nap.min() if len(dists_nap) > 0 else np.nan

    # CHI–PAR minimum distance
    pairs_par, dists_par = capped_distance(
        sel_chi.positions, sel_par.positions,
        max_cutoff=100.0, box=ts.dimensions
    )
    dmin_par = dists_par.min() if len(dists_par) > 0 else np.nan

    records.append([ts.frame, ts.time, dmin_nap, dmin_par])

df = pd.DataFrame(records, columns=["frame", "time_ps", "dmin_chi_nap", "dmin_chi_par"])
df.to_csv("min_distances.csv", index=False)

print("Output: min_distances.csv")
print(f"Mean min CHI–NAP distance: {df['dmin_chi_nap'].mean():.2f} Å")
print(f"Mean min CHI–PAR distance: {df['dmin_chi_par'].mean():.2f} Å