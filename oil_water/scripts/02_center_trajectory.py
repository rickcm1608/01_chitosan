"""
Center and wrap the production trajectory so that chitosan stays in the
center of the box, removing PBC artifacts before analysis.

Requires: MDAnalysis
Input:  complex.prmtop, traj_50ns.dcd
Output: traj_centered_50ns.dcd
"""

import MDAnalysis as mda
from MDAnalysis.transformations import wrap, center_in_box

# ===============================
# Load
# ===============================
u = mda.Universe("complex.prmtop", "traj_50ns.dcd")

# Select chitosan as the centering group
chi = u.select_atoms("resname CHI")
all_atoms = u.atoms

# ===============================
# Apply on-the-fly transformations
# ===============================
transforms = [
    center_in_box(chi, wrap=True),
    wrap(all_atoms)
]
u.trajectory.add_transformations(*transforms)

# ===============================
# Write centered trajectory
# ===============================
print("Writing centered trajectory...")
with mda.Writer("traj_centered_50ns.dcd", n_atoms=u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)

print("Output: traj_centered_50ns.dcd")
