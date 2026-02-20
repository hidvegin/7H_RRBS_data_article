#!/usr/bin/env python3
"""
da04_chromosome_coverage.py
---------------------------
Visualises the distribution of methylation probes across chromosomes
for each context and sample.

Outputs:
  chromosome_coverage.png/.svg   - grouped bar chart: probes per chromosome × sample
  chromosome_coverage_summary.csv - already produced by da01, but this adds % of total
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..')
OUT_DIR    = SCRIPT_DIR

CONTEXTS = ['CpG', 'CHG', 'CHH']
SAMPLES  = ['Mv9', 'Igri', '7Hadd']

SAMPLE_COLORS = {
    'Mv9':   '#2166AC',
    'Igri':  '#D6604D',
    '7Hadd': '#4DAC26',
}

# Ordered chromosome list for plotting
WHEAT_CHROMS  = [f'{i}{s}' for i in range(1, 8) for s in ['A', 'B', 'D']]
BARLEY_CHROMS = [f'{i}H' for i in range(1, 8)]
CHROM_ORDER   = WHEAT_CHROMS + BARLEY_CHROMS

print("=== da04: Chromosome Coverage ===")

fig, axes = plt.subplots(3, 1, figsize=(20, 14), sharex=False)
fig.suptitle('Methylation Probe Coverage per Chromosome', fontsize=15, fontweight='bold')

all_rows = []

for ax, ctx in zip(axes, CONTEXTS):
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    print(f"Loading {ctx}...")
    df = pd.read_csv(fpath, sep='\t',
                     usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                     low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)

    # Count non-NaN sites per chromosome per sample
    records = []
    for chrom, grp in df.groupby('Chromosome'):
        row = {'Context': ctx, 'Chromosome': chrom}
        for s in SAMPLES:
            row[s] = int(pd.to_numeric(grp[s], errors='coerce').notna().sum())
        row['Total_probes'] = len(grp)
        records.append(row)
        all_rows.append(row)

    cov_df = pd.DataFrame(records)
    # Keep only chromosomes we know about; sort them
    cov_df = cov_df[cov_df['Chromosome'].isin(CHROM_ORDER)].copy()
    cov_df['_order'] = cov_df['Chromosome'].map({c: i for i, c in enumerate(CHROM_ORDER)})
    cov_df = cov_df.sort_values('_order').drop(columns='_order')

    chroms = cov_df['Chromosome'].tolist()
    x = np.arange(len(chroms))
    width = 0.25

    for j, sample in enumerate(SAMPLES):
        ax.bar(x + j * width, cov_df[sample], width,
               label=sample, color=SAMPLE_COLORS[sample], alpha=0.85, edgecolor='white')

    # Vertical line separating wheat from barley
    wheat_count = sum(c in WHEAT_CHROMS for c in chroms)
    ax.axvline(wheat_count - 0.5 + 1.5 * width, color='black', ls='--', lw=1.2, alpha=0.6)
    ax.text(wheat_count - 0.2 + 1.5 * width, ax.get_ylim()[1] * 0.02,
            '← Wheat   Barley →', fontsize=8, ha='center', color='black', alpha=0.7)

    ax.set_title(f'{ctx} context', fontsize=12)
    ax.set_ylabel('Sites with data', fontsize=10)
    ax.set_xticks(x + width)
    ax.set_xticklabels(chroms, rotation=60, ha='right', fontsize=8)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'chromosome_coverage.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")
plt.close()

# Save extended summary CSV
summary_df = pd.DataFrame(all_rows)
summary_df.to_csv(os.path.join(OUT_DIR, 'chromosome_coverage_detail.csv'), index=False)
print(f"Saved chromosome_coverage_detail.csv ({len(summary_df)} rows)")

print("\n✓ da04 complete.")
