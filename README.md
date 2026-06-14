# Chitosan Oil-Water MD Simulation

**Paper:** "Chitosan Eco-Friendly Approach to Oil Spill Cleanup: A Combined 2D TD-NMR Relaxation and Computational Modeling Study"
**Journal:** Journal of Polymers and the Environment (2025)

## Overview

Molecular dynamics study of chitosan as an oil-spill remediation agent, comparing two systems:

| Folder | System | Water model |
|--------|--------|-------------|
| `water_only/` | Chitosan in TIP4PEW water (267 CHI atoms, 49 Na+, 49 Cl-) | TIP4PEW |
| `oil_water/` | Chitosan + 20 naphthalene + 20 paraffin in TIP4PEW water (60 Na+, 60 Cl-) | TIP4PEW |

Production run: **50 ns** (25,000,000 steps x 2 fs, 1000 frames saved every 50 ps)
Software: OpenMM (simtk API) + AMBER for minimization/equilibration

---

## Repository Structure

```
chitosan-oil-water-md_github/
+-- parametrization/
|   +-- chitosan/          # GAFF2: chitosan.frcmod, .mol2, .pdb, sqm files
|   +-- naphthalene/       # GAFF2: naphthalene.frcmod, .mol2, .pdb
|   +-- paraffin/          # GAFF2: paraffin.frcmod, .mol2, .pdb
+-- water_only/
|   +-- topology/          # complex.prmtop, complex.inpcrd, tleap.in
|   +-- amber_protocol/    # min.in, heat.in, equil1-4.in
|   +-- scripts/           # Analysis scripts (see below)
|   +-- data/
|       +-- equilibration/ # .out files and parsed CSVs (energy, temp, density)
|       +-- rdf_blocks/    # Per-block CHI-WAT RDF (10 blocks x 5 ns)
|       +-- rdf_averaged/  # Equilibrated block average (blocks 5-10)
+-- oil_water/
    +-- topology/          # complex.prmtop, complex.inpcrd, tleap.in
    +-- amber_protocol/    # min.in, heat.in, equil1-4.in
    +-- scripts/           # Analysis scripts (see below)
    +-- data/
        +-- equilibration/ # .out files and parsed CSVs
        +-- rdf_chi_wat/   # Per-block CHI-WAT RDF (25 blocks x 2 ns)
        +-- rdf_chi_nap/   # Per-block CHI-NAP RDF
        +-- rdf_chi_par/   # Per-block CHI-PAR RDF
```

---

## Scripts

### water_only/scripts/

| Script | Description | Output |
|--------|-------------|--------|
| `01_production_md.py` | 50 ns NPT production (OpenMM, CUDA) | `traj_50ns.dcd` |
| `02_water_coordination.py` | Water O atoms within 3.5 A of CHI per frame | `water_coordination.csv` |
| `03_hbonds_chi_water.py` | H-bond count CHI-WAT per frame (distance <= 3.5 A) | `hbonds_chi_wat.csv` |
| `04_rdf_block_convergence.py` | CHI-WAT RDF over 10 blocks (convergence check) | `rdf_blocks/` |
| `05_rdf_equilibrated_average.py` | Block average of equilibrated RDF (blocks 5-10) | `rdf_averaged/chi_wat_rdf_avg.txt` |
| `06_radius_of_gyration.py` | Rg of all atoms per frame | `radius_of_gyration.csv` |

### oil_water/scripts/

| Script | Description | Output |
|--------|-------------|--------|
| `01_production_md.py` | 50 ns NPT production (OpenMM, CUDA) | `traj_50ns.dcd` |
| `02_center_trajectory.py` | Center CHI in box, wrap PBC | `traj_centered_50ns.dcd` |
| `03_parse_equilibration.py` | Parse AMBER .out files for QC | `minimization.csv`, `heating.csv`, `equil1_nvt.csv`, `equil2_npt.csv` |
| `04_rdf_block_convergence.py` | CHI-WAT/NAP/PAR RDF over 25 blocks | `rdf_chi_wat/`, `rdf_chi_nap/`, `rdf_chi_par/` |
| `05_rdf_equilibrated_average.py` | Block average of last 4 blocks per pair | `rdf_averaged/<pair>_rdf_avg.csv` |
| `06_hbonds_chi_water.py` | H-bond count CHI-WAT per frame | `hbonds_chi_wat.csv` |
| `07_water_coordination.py` | Water coordination number around CHI | `water_coordination.csv` |
| `08_com_distances.py` | COM distances CHI-NAP and CHI-PAR | `com_distances.csv` |
| `09_min_distances.py` | Minimum pairwise distances CHI-NAP and CHI-PAR | `min_distances.csv` |
| `10_radius_of_gyration.py` | Rg of all atoms per frame | `radius_of_gyration.csv` |

---

## System Preparation

Force field assignments:
- Chitosan: GAFF2 + custom `chitosan.frcmod` (antechamber/sqm)
- Naphthalene: GAFF2 + `naphthalene.frcmod`
- Paraffin: GAFF2 + `paraffin.frcmod`
- Water: TIP4PEW (`leaprc.water.tip4pew`, TIP4PEWBOX)
- Ions: Na+, Cl- (neutralization + 150 mM NaCl)

AMBER equilibration protocol (input files in each system's `amber_protocol/`):
1. Energy minimization (`min.in`)
2. NVT heating 0->300 K (`heat.in`)
3. NVT equilibration at 300 K (`equil1.in`, `equil2.in`)
4. NPT equilibration at 300 K, 1 atm (`equil3.in`, `equil4.in`)

---

## MD Protocol (OpenMM)

```
Ensemble:         NPT (Langevin + Monte Carlo Barostat)
Temperature:      300 K  (gamma = 1 ps-1)
Pressure:         1 atm  (barostat every 25 steps)
Time step:        2 fs
Constraints:      HBonds
Nonbonded:        PME, cutoff 1.0 nm
Production:       50 ns (25 M steps)
Save frequency:   every 50 ps -> 1000 frames
GPU:              CUDA
```

---

## Dependencies

```
python >= 3.8
openmm (simtk API)
mdanalysis >= 2.0
numpy
pandas
matplotlib
```

---

## Key Results

- **Rg (water-only):** ~35.5 A -- chitosan adopts a compact conformation in pure water
- **Rg (oil-water):** ~41.9 A -- chitosan adopts a more extended conformation upon interaction with oil molecules
- **RDF CHI-WAT:** first coordination shell at ~3.7 A; reduced in oil-water system
- **RDF CHI-NAP/PAR:** direct adsorption contacts confirmed at short range

Data used in paper figures: `rdf_averaged/chi_wat_rdf_avg.txt` (water_only) and `rdf_chi_wat/` blocks 21-25 average (oil_water); `radius_of_gyration.csv` from both systems.

## Citation

If you use this code or data, please cite:

> Kock, F.V.C.; Valdiviezo, J.; Cirilo, E.; Sifuentes, J.; Nakamatsu, J.; Souza, A.; Castro, E.; Barbosa, L.; Colnago, L. "Chitosan Eco-Friendly Approach to Oil Spill Cleanup: A Combined 2D TD-NMR Relaxation and Computational Modeling Study." *Journal of Polymers and the Environment* (2025). https://doi.org/10.1007/s10924-025-03722-1
