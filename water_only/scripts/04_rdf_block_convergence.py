import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import numpy as np
import matplotlib.pyplot as plt
import os

# -----------------------------
# Input files
# -----------------------------
top  = "complex.prmtop"
traj = "traj_centered_50ns.dcd"

u = mda.Universe(top, traj)

# -----------------------------
# Atom selections
# -----------------------------
sel_chi = u.select_atoms("resname CHI")                    # chitosan
sel_wat = u.select_atoms("resname WAT and name O")         # water oxygens

# -----------------------------
# RDF parameters
# -----------------------------
r_max = 40.0   # maximum r in Angstroms
nbins = 200

# -----------------------------
# Block decomposition for convergence check
# -----------------------------
n_frames  = len(u.trajectory)
n_blocks  = 10
block_size = n_frames // n_blocks

outdir = "rdf_blocks"
os.makedirs(outdir, exist_ok=True)

plt.figure(figsize=(8, 6))

for i in range(n_blocks):
    start = i * block_size
    end   = (i + 1) * block_size

    rdf_block = InterRDF(sel_chi, sel_wat, range=(0.0, r_max), nbins=nbins)
    rdf_block.run(start=start, stop=end)

    r = rdf_block.results.bins
    g = rdf_block.results.rdf

    plt.plot(r, g, label=f"Block {i+1}")

    out_file = os.path.join(outdir, f"chi_wat_block_{i+1:02d}.txt")
    np.savetxt(out_file, np.column_stack([r, g]), header="r[Å] g(r)")

plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title("CHI–WAT RDF block convergence (10 blocks × 5 ns)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "chi_wat_convergence.png"), dpi=300)
plt.show()

print(f"Block RDF