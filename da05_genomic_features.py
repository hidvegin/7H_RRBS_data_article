#!/usr/bin/env python3
"""
da05_genomic_features.py
------------------------
Analyses the distribution of differentially methylated probes across
genomic features (promoter / gene body / downstream) broken down by
methylation context and genome-of-origin (Wheat vs. Barley 7H).

Input:  ../filtered_significant_methylation.csv  (SeqMonk-significant probes)
Outputs:
  genomic_features_summary.csv      - counts & proportions per feature × context
  genomic_features_barplot.png/.svg - stacked bar chart
  hyper_hypo_breakdown.csv          - hyper- vs. hypomethylation per context × region
  hyper_hypo_barplot.png/.svg       - stacked bar chart for hyper/hypo
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

INPUT_FILE = os.path.join(DATA_DIR, 'filtered_significant_methylation.csv')

print("=== da05: Genomic Feature Distribution ===")
df = pd.read_csv(INPUT_FILE, low_memory=False)
print(f"Loaded {len(df):,} significant probes.")

# ── Clean up region and status labels ────────────────────────────────────────
REGION_MAP = {
    'overlapping': 'Gene body',
    'upstream':    'Promoter',
    'downstream':  'Downstream',
}
df['Region_label'] = df['Region'].map(REGION_MAP).fillna(df['Region'])

# Comparison: split into Wheat (Mv9 vs 7Hadd) and Barley (Igri vs 7Hadd)
df['Genome'] = df['Comparison'].apply(
    lambda x: 'Barley 7H' if 'Barley' in str(x) or '7H' in str(x) else 'Wheat (Mv9)'
)

CONTEXTS = ['CpG', 'CHG', 'CHH']
REGIONS  = ['Promoter', 'Gene body', 'Downstream']
GENOMES  = ['Wheat (Mv9)', 'Barley 7H']

# ── 1. Feature × Context counts ──────────────────────────────────────────────

summary_rows = []
for genome in GENOMES:
    sub = df[df['Genome'] == genome]
    for ctx in CONTEXTS:
        ctx_sub = sub[sub['Context'] == ctx]
        total = len(ctx_sub)
        for reg in REGIONS:
            count = (ctx_sub['Region_label'] == reg).sum()
            summary_rows.append({
                'Genome': genome,
                'Context': ctx,
                'Region': reg,
                'Count': count,
                'Pct_of_context': round(count / total * 100, 2) if total > 0 else 0,
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT_DIR, 'genomic_features_summary.csv'), index=False)
print(summary_df.to_string(index=False))

# ── 2. Stacked bar chart: feature distribution ─────────────────────────────

REGION_COLORS = {
    'Promoter':   '#E7298A',
    'Gene body':  '#1B7837',
    'Downstream': '#7570B3',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
fig.suptitle('Genomic Feature Distribution of Differentially Methylated Probes',
             fontsize=13, fontweight='bold')

for ax, genome in zip(axes, GENOMES):
    sub = summary_df[summary_df['Genome'] == genome].copy()
    pivot = sub.pivot(index='Context', columns='Region', values='Count').reindex(
        index=CONTEXTS, columns=REGIONS, fill_value=0
    )

    bottoms = np.zeros(len(CONTEXTS))
    x = np.arange(len(CONTEXTS))
    for reg in REGIONS:
        vals = pivot[reg].values
        ax.bar(x, vals, bottom=bottoms, label=reg,
               color=REGION_COLORS[reg], alpha=0.88, edgecolor='white')
        # Label inside bar if large enough
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v > 50:
                ax.text(xi, b + v / 2, str(int(v)),
                        ha='center', va='center', fontsize=9, color='white', fontweight='bold')
        bottoms += vals

    ax.set_title(f'{genome}', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(CONTEXTS, fontsize=11)
    ax.set_ylabel('Number of probes', fontsize=10)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    ax.legend(title='Region', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'genomic_features_barplot.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")
plt.close()

# ── 3. Hyper / Hypo breakdown ─────────────────────────────────────────────────

hh_rows = []
for genome in GENOMES:
    sub = df[df['Genome'] == genome]
    # Which status column to use
    status_col = 'Igri_7Hadd_Status' if genome == 'Barley 7H' else 'Mv9_7Hadd_Status'
    for ctx in CONTEXTS:
        ctx_sub = sub[sub['Context'] == ctx]
        for reg in REGIONS:
            reg_sub = ctx_sub[ctx_sub['Region_label'] == reg]
            hyper = (reg_sub[status_col].str.lower() == 'hyper').sum()
            hypo  = (reg_sub[status_col].str.lower() == 'hypo').sum()
            hh_rows.append({
                'Genome': genome, 'Context': ctx, 'Region': reg,
                'Hyper': int(hyper), 'Hypo': int(hypo), 'Total': int(hyper + hypo),
            })

hh_df = pd.DataFrame(hh_rows)
hh_df.to_csv(os.path.join(OUT_DIR, 'hyper_hypo_breakdown.csv'), index=False)

# Plot hyper/hypo stacked bars per genome × context
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=False)
fig.suptitle('Hyper- vs. Hypomethylation by Genomic Region', fontsize=13, fontweight='bold')

for row_idx, genome in enumerate(GENOMES):
    for col_idx, ctx in enumerate(CONTEXTS):
        ax = axes[row_idx][col_idx]
        sub = hh_df[(hh_df['Genome'] == genome) & (hh_df['Context'] == ctx)]
        x = np.arange(len(REGIONS))
        hyper_vals = sub.set_index('Region').reindex(REGIONS)['Hyper'].fillna(0).values
        hypo_vals  = sub.set_index('Region').reindex(REGIONS)['Hypo'].fillna(0).values

        ax.bar(x, hyper_vals, label='Hyper', color='#D73027', alpha=0.85, edgecolor='white')
        ax.bar(x, hypo_vals, bottom=hyper_vals, label='Hypo', color='#4575B4', alpha=0.85, edgecolor='white')

        for xi, (h, l) in enumerate(zip(hyper_vals, hypo_vals)):
            if h > 5:
                ax.text(xi, h / 2, str(int(h)), ha='center', va='center',
                        fontsize=8, color='white', fontweight='bold')
            if l > 5:
                ax.text(xi, h + l / 2, str(int(l)), ha='center', va='center',
                        fontsize=8, color='white', fontweight='bold')

        ax.set_title(f'{genome} – {ctx}', fontsize=10)
        ax.set_xticks(x)
        ax.set_xticklabels([r[:4] for r in REGIONS], fontsize=9)
        ax.set_ylabel('Probes', fontsize=9)
        if col_idx == 0 and row_idx == 0:
            ax.legend(fontsize=9)
        ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'hyper_hypo_barplot.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")
plt.close()

print("\n✓ da05 complete.")
