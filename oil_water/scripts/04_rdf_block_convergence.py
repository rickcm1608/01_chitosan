"""
Block-decomposed RDF analysis for convergence checking.

Computes CHI–WAT, CHI–NAP, and CHI–PAR radial distribution functions
over 25 equal-length blocks (each ~2 ns for a 50 ns trajectory).
Each pair is saved to its own subfolder.

Input:  complex.prmtop, traj_centered_50ns.dcd
Output: rdf_chi_wat/block_<nn>.txt  + rdf_chi_wat/convergence.png
        rdf_chi_nap/block_<nn>.txt  + rdf_chi_nap/convergence.png
        rdf_chi_par/block_<nn>.txt  + rdf_chi_par/convergence.png
"""

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import matplotlib.pyplot as plt
import os

# -----------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------
u = mda.Universe("complex.prmtop", "traj_centered_50ns.dcd")

sel_chi = u.select_atoms("resname CHI")
sel_wat = u.select_atoms("resname WAT and name O")
sel_nap = u.select_atoms("resname NAP")
sel_par = u.select_atoms("resname PAR")

pairs = {
    "chi_wat": (sel_chi, sel_wat),
    "chi_nap": (sel_chi, sel_nap),
    "chi_par": (sel_chi, sel_par),
}

# -----------------------------------------------------------------------
# Block parameters
# -----------------------------------------------------------------------
r_max      = 40.0
nbins      = 200
n_blocks   = 25
n_frames   = len(u.trajectory)
block_size = n_frames // n_blocks

# -----------------------------------------------------------------------
# Compute RDF for each pair × block — one output folder per pair
# -----------------------------------------------------------------------
for pair_name, (sel_a, sel_b) in pairs.items():
    outdir = f"rdf_{pair_name}"
    os.makedirs(outdir, exist_ok=True)

    plt.figure(figsize=(8, 5))

    for i in range(n_blocks):
        start = i * block_size
        end   = (i + 1) * block_size

        rdf = InterRDF(sel_a, sel_b, range=(0.0, r_max), nbins=nbins)
        rdf.run(start=start, stop=end)

        r = rdf.results.bins
        g = rdf.results.rdf

        plt.plot(r, g, alpha=0.6, label=f"Block {i+1}")

        out_file = os.path.join(outdir, f"block_{i+1:02d}.txt")
        np.savetxt(out_file, np.column_stack([r, g]), header="r[Å] g(r)")

    label = pair_name.replace("_", "–").upper()
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"{label} RDF — 25 blocks × 2 ns")
    plt.tight_layout()
    plt.sa