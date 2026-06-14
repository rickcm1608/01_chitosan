"""
Radius of gyration over the production trajectory.

Input:  complex.prmtop, traj_centered_50ns.dcd
Output: radius_of_gyration.csv  (frame, rg_ang)
"""

import MDAnalysis as mda
import pandas as pd

# ===============================
# Load
# ===============================
u = mda.Universe("complex.prmtop", "traj_centered_50ns.dcd")

rg_list = []

print("Computing radius of gyration per frame...")
for ts in u.trajectory:
    rg = u.atoms.radius_of_gyration()
    rg_list.append(rg)

df = pd.DataFrame({"frame": range(len(rg_list)), "rg_ang": rg_list})
df.to_csv("radius_of_gyration.csv", index=False)

print("Output: radius_of_gyration.csv")
print(f"Mean Rg: {df[