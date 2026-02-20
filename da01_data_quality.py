#!/usr/bin/env python3
"""
da01_data_quality.py
--------------------
Data quality and coverage assessment for Scientific Data submission.

Outputs (all written to this script's directory):
  data_quality_summary.csv       - site counts and methylation statistics per context × sample
  data_completeness.csv          - sites present in 1/2/3 samples per context
  bisulfite_conversion_estimate.csv - mean CHH methylation as conversion proxy
"""

import os
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..')
OUT_DIR    = SCRIPT_DIR

CONTEXTS = ['CpG', 'CHG', 'CHH']
SAMPLES  = ['Mv9', 'Igri', '7Hadd']

WHEAT_CHROMS = [f'{i}{s}' for i in range(1, 8) for s in ['A', 'B', 'D']]
BARLEY_CHROMS = [f'{i}H' for i in range(1, 8)]

# ── 1. Global quality summary ────────────────────────────────────────────────

print("=== da01: Data Quality Summary ===")
quality_rows = []
completeness_rows = []

for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    print(f"\nLoading {ctx}...")
    df = pd.read_csv(fpath, sep='\t',
                     usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                     low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)

    total_sites = len(df)

    # Per-sample stats
    for sample in SAMPLES:
        vals = pd.to_numeric(df[sample], errors='coerce')
        n_total    = total_sites
        n_with_data = vals.notna().sum()
        quality_rows.append({
            'Context': ctx,
            'Sample': sample,
            'Total_sites_in_matrix': n_total,
            'Sites_with_data': int(n_with_data),
            'Coverage_pct': round(n_with_data / n_total * 100, 2),
            'Mean_methylation': round(vals.dropna().mean(), 4),
            'Median_methylation': round(vals.dropna().median(), 4),
            'Std_methylation': round(vals.dropna().std(), 4),
            'Min_methylation': round(vals.dropna().min(), 4),
            'Max_methylation': round(vals.dropna().max(), 4),
        })

    # Completeness: how many sites have data in 1 / 2 / all 3 samples
    has_data = df[SAMPLES].notna()
    n_samples_with_data = has_data.sum(axis=1)
    completeness_rows.append({
        'Context': ctx,
        'Total_sites': total_sites,
        'Sites_in_all_3_samples': int((n_samples_with_data == 3).sum()),
        'Sites_in_exactly_2_samples': int((n_samples_with_data == 2).sum()),
        'Sites_in_exactly_1_sample': int((n_samples_with_data == 1).sum()),
        'Complete_case_pct': round((n_samples_with_data == 3).sum() / total_sites * 100, 2),
    })

quality_df = pd.DataFrame(quality_rows)
completeness_df = pd.DataFrame(completeness_rows)

quality_df.to_csv(os.path.join(OUT_DIR, 'data_quality_summary.csv'), index=False)
completeness_df.to_csv(os.path.join(OUT_DIR, 'data_completeness.csv'), index=False)

print("\n── Quality Summary ──")
print(quality_df.to_string(index=False))
print("\n── Completeness ──")
print(completeness_df.to_string(index=False))

# ── 2. Bisulfite conversion efficiency estimate from CHH context ─────────────

print("\n\n=== Bisulfite Conversion Estimate (from CHH) ===")
chh_path = os.path.join(DATA_DIR, 'clean_CHH_all.txt')
chh = pd.read_csv(chh_path, sep='\t',
                  usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                  low_memory=False)
chh['Chromosome'] = chh['Chromosome'].astype(str)

# Use only wheat chromosomes for the conversion estimate
# (barley CHH biology differs; wheat non-CG meth in genes should be minimal)
chh_wheat = chh[chh['Chromosome'].isin(WHEAT_CHROMS)]

conv_rows = []
for sample in SAMPLES:
    vals = pd.to_numeric(chh_wheat[sample], errors='coerce').dropna()
    # Conversion efficiency = 100 - mean_CHH (approximation)
    mean_chh = vals.mean()
    conv_est = 100 - mean_chh
    conv_rows.append({
        'Sample': sample,
        'Mean_CHH_methylation_pct': round(mean_chh, 4),
        'Estimated_conversion_efficiency_pct': round(conv_est, 4),
        'Note': 'CHH on wheat chromosomes only; lower mean CHH = better conversion'
    })

conv_df = pd.DataFrame(conv_rows)
conv_df.to_csv(os.path.join(OUT_DIR, 'bisulfite_conversion_estimate.csv'), index=False)
print(conv_df.to_string(index=False))

# ── 3. Per-chromosome site counts ────────────────────────────────────────────

print("\n\n=== Per-Chromosome Site Counts ===")
chrom_rows = []
for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    df = pd.read_csv(fpath, sep='\t',
                     usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                     low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)
    for chrom, grp in df.groupby('Chromosome'):
        row = {'Context': ctx, 'Chromosome': chrom}
        for s in SAMPLES:
            row[f'{s}_sites'] = int(grp[s].notna().sum())
        row['Total_probes'] = len(grp)
        chrom_rows.append(row)

chrom_df = pd.DataFrame(chrom_rows)
chrom_df.to_csv(os.path.join(OUT_DIR, 'per_chromosome_site_counts.csv'), index=False)
print(f"Saved {len(chrom_df)} chromosome × context rows.")

print("\n✓ da01 complete. Outputs written to:", OUT_DIR)
