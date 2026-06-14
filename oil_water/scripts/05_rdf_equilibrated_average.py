"""
Equilibrated RDF block average.

Uses the last 4 of 25 blocks (blocks 22–25, ~8 ns of equilibrated data)
to compute final CHI–WAT, CHI–NAP, and CHI–PAR RDFs with error estimates.

Reads block files produced by 04_rdf_block_convergence.py.

Input:  rdf_chi_wat/block_<nn>.txt, rdf_chi_nap/..., rdf_chi_par/...
Output: rdf_averaged/chi_wat_rdf_avg.csv  +  rdf_averaged/chi_wat_rdf_avg.png
        rdf_averaged/chi_nap_rdf_avg.csv  +  rdf_averaged/chi_nap_rdf_avg.png
        rdf_averaged/chi_par_rdf_avg.csv  +  rdf_averaged/chi_par_rdf_avg.png
"""

import numpy as np
import matplotlib.pyplot as plt
import os

pairs    = ["chi_wat", "chi_nap", "chi_par"]
n_blocks = 25
last_k   = 4          # average last 4 blocks (22–25)
outdir   = "rdf_averaged"

os.makedirs(outdir, exist_ok=True)

for pair in pairs:
    indir    = f"rdf_{pair}"
    rdf_list = []
    r_values = None

    for i in range(n_blocks - last_k, n_blocks):   # blocks 22–25 (0-indexed)
        fname = os.path.join(indir, f"block_{i+1:02d}.txt")
        data  = np.loadtxt(fname)
        if r_values is None:
            r_values = data[:, 0]
        rdf_list.append(data[:, 1])

    rdf_array = np.array(rdf_list)
    g_mean    = rdf_array.mean(axis=0)
    g_std     = rdf_array.std(axis=0)

    # Save CSV
    out_csv = os.path.join(outdir, f"{pair}_rdf_avg.csv")
    np.savetxt(
        out_csv,
        np.column_stack([r_values, g_mean, g_std]),
        header="r_ang,g_mean,g_std",
        delimiter=","
    )

    # Figure with error band
    label = pair.replace("_", "–").upper()
    plt.figure(figsize=(8, 5))
    plt.plot(r_values, g_mean, label=f"{label} (blocks 22–25)")
    plt.fill_between(r_values, g_mean - g_std, g_mean + g_std, alpha=0.3)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"{label} RDF — block average (equilibrated)")
    plt.legend()
    plt.tight_layout()
 