"""
Parse AMBER equilibration output files and extract thermodynamic properties
for quality control.

Input files (from AMBER minimization + equilibration protocol):
  min.out    → energy minimization
  heat.out   → NVT heating (0 → 300 K)
  equil1.out → NVT equilibration at 300 K
  equil2.out → NPT equilibration at 300 K, 1 atm

Output CSV files (one per stage):
  minimization.csv, heating.csv, equil1_nvt.csv, equil2_npt.csv
"""

import re
import pandas as pd


# -----------------------------------------------------------------------
# Helper parsers
# -----------------------------------------------------------------------

def parse_min(filename):
    """Parse AMBER minimization output → Step, Energy, RMS, GMAX."""
    records = []
    with open(filename) as f:
        text = f.read()
    pattern = re.compile(
        r'NSTEP\s*=\s*(\d+)\s+ENERGY\s*=\s*([\d.E+\-]+)'
        r'.*?RMS\s*=\s*([\d.E+\-]+)\s+GMAX\s*=\s*([\d.E+\-]+)'
        r'.*?Atom\s*=\s*(\d+)',
        re.DOTALL
    )
    for m in pattern.finditer(text):
        records.append({
            'step':        int(m.group(1)),
            'energy':      float(m.group(2)),
            'rms':         float(m.group(3)),
            'gmax':        float(m.group(4)),
            'atom_number': int(m.group(5)),
        })
    return pd.DataFrame(records)


def parse_nvt(filename):
    """Parse AMBER NVT output (heating or equil1) → temp_k, etot."""
    records = []
    with open(filename) as f:
        text = f.read()
    pattern = re.compile(
        r'NSTEP\s*=\s*\d+.*?'
        r'TEMP\(K\)\s*=\s*([\d.E+\-]+).*?'
        r'Etot\s*=\s*([\d.E+\-]+)',
        re.DOTALL
    )
    for m in pattern.finditer(text):
        records.append({
            'temp_k': float(m.group(1)),
            'etot':   float(m.group(2)),
        })
    return pd.DataFrame(records)


def parse_npt(filename):
    """Parse AMBER NPT output (equil2) → nstep, time_ps, temp_k, press, energies, volume, density."""
    records = []
    with open(filename) as f:
        text = f.read()
    pattern = re.compile(
        r'NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([\d.E+\-]+).*?'
        r'TEMP\(K\)\s*=\s*([\d.E+\-]+)\s+PRESS\s*=\s*([\d.E+\-]+).*?'
        r'Etot\s*=\s*([\d.E+\-]+)\s+EKtot\s*=\s*([\d.E+\-]+)\s+EPtot\s*=\s*([\d.E+\-]+).*?'
        r'VOLUME\s*=\s*([\d.E+\-]+).*?'
        r'Density\s*=\s*([\d.E+\-]+)',
        re.DOTALL
    )
    for m in pattern.finditer(text):
        records.append({
            'nstep':   int(m.group(1)),
            'time_ps': float(m.group(2)),
            'temp_k':  float(m.group(3)),
            'press':   float(m.group(4)),
            'etot':    float(m.group(5)),
            'ektot':   float(m.group(6)),
            'eptot':   float(m.group(7)),
            'volume':  float(m.group(8)),
            'density': float(m.group(9)),
        })
    return pd.DataFrame(records)


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

df_min   = parse_min('min.out')
df_heat  = parse_nvt('heat.out')
df_eq1   = parse_nvt('equil1.out')
df_eq2   = parse_npt('equil2.out')

df_min.to_csv('minimization.csv',  index=False)
df_heat.to_csv('heating.csv',      index=False)
df_eq1.to_csv('equil1_nvt.csv',    index=False)
df_eq2.to_csv('equil2_npt.csv',