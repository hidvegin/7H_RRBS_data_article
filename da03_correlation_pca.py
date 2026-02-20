#!/usr/bin/env python3
"""
da03_correlation_pca.py
-----------------------
Pairwise sample correlations and chromosome-level PCA for Scientific Data.

Correlation strategy:
  - Mv9 vs 7Hadd : wheat chromosomes (1A-7D), sites with data in BOTH samples
  - Igri vs 7Hadd: 7H chromosome only, sites with data in BOTH samples
  - Mv9 vs Igri  : 7H chromosome (shared reference region), sites with data in BOTH

Note: very few "all-3" complete cases exist (~400–500 per context) because
the three samples represent different genome compositions
(wheat-only / barley-only / wheat+7H hybrid), so pairwise analysis is correct.

PCA: chromosome-level (mean methylation per chromosome per sample per context).
This avoids the sparsity problem while providing a meaningful overview.

Outputs:
  correlation_summary.csv          - Pearson & Spearman r per context × sample pair
  correlation_heatmaps.png/.svg    - 3×2 heatmap grid
  pca_chromosome_level.png/.svg    - chromosome-level PCA
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from itertools import combinations

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..')
OUT_DIR    = SCRIPT_DIR

CONTEXTS = ['CpG', 'CHG', 'CHH']
SAMPLES  = ['Mv9', 'Igri', '7Hadd']

WHEAT_CHROMS  = [f'{i}{s}' for i in range(1, 8) for s in ['A', 'B', 'D']]
BARLEY_7H     = ['7H']

SAMPLE_COLORS = {
    'Mv9':   '#2166AC',
    'Igri':  '#D6604D',
    '7Hadd': '#4DAC26',
}

CONTEXT_COLORS = {'CpG': '#1A9641', 'CHG': '#FDAE61', 'CHH': '#D7191C'}

print("=== da03: Correlation & PCA ===")

# ── 1. Pairwise correlations (biologically meaningful pairs) ──────────────────

# Pair definitions: (sample_A, sample_B, chromosomes_to_use, label)
PAIRS = [
    ('Mv9',  '7Hadd', WHEAT_CHROMS, 'Mv9 vs 7Hadd (wheat chr)'),
    ('Igri', '7Hadd', BARLEY_7H,    'Igri vs 7Hadd (7H chr)'),
    ('Mv9',  'Igri',  WHEAT_CHROMS, 'Mv9 vs Igri (wheat chr)'),
]

corr_rows = []
# For heatmap: store per-context full 3×3 matrices (using pairwise-available data)
pearson_mats  = {ctx: pd.DataFrame(np.nan, index=SAMPLES, columns=SAMPLES) for ctx in CONTEXTS}
spearman_mats = {ctx: pd.DataFrame(np.nan, index=SAMPLES, columns=SAMPLES) for ctx in CONTEXTS}
for ctx in CONTEXTS:
    for s in SAMPLES:
        pearson_mats[ctx].loc[s, s]  = 1.0
        spearman_mats[ctx].loc[s, s] = 1.0

for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    print(f"\nLoading {ctx}...")
    df = pd.read_csv(fpath, sep='\t',
                     usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                     low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)
    for col in SAMPLES:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    for s1, s2, chroms, label in PAIRS:
        sub = df[df['Chromosome'].isin(chroms)][[s1, s2]].dropna()
        n = len(sub)
        if n < 10:
            print(f"  SKIP {label} ({ctx}): only {n} shared sites")
            continue
        r_p, pval_p = stats.pearsonr(sub[s1], sub[s2])
        r_s, pval_s = stats.spearmanr(sub[s1], sub[s2])
        corr_rows.append({
            'Context': ctx,
            'Sample_A': s1, 'Sample_B': s2,
            'Chromosomes': 'wheat' if chroms == WHEAT_CHROMS else '7H',
            'N_shared_sites': n,
            'Pearson_r': round(r_p, 6),
            'Pearson_p': f'{pval_p:.2e}',
            'Spearman_r': round(r_s, 6),
            'Spearman_p': f'{pval_s:.2e}',
            'Label': label,
        })
        pearson_mats[ctx].loc[s1, s2]  = round(r_p, 4)
        pearson_mats[ctx].loc[s2, s1]  = round(r_p, 4)
        spearman_mats[ctx].loc[s1, s2] = round(r_s, 4)
        spearman_mats[ctx].loc[s2, s1] = round(r_s, 4)
        print(f"  {label}: n={n:,}  Pearson r={r_p:.4f}  Spearman r={r_s:.4f}")

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(os.path.join(OUT_DIR, 'correlation_summary.csv'), index=False)
print("\n── Correlation Summary saved ──")

# ── 2. Correlation heatmaps ───────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
fig.suptitle(
    'Pairwise Sample Correlations by Methylation Context\n'
    '(Mv9↔7Hadd & Mv9↔Igri: wheat chr; Igri↔7Hadd: 7H chr)',
    fontsize=12, fontweight='bold'
)

titles = [('Pearson r', pearson_mats), ('Spearman r', spearman_mats)]
for row_idx, (method, mats) in enumerate(titles):
    for col_idx, ctx in enumerate(CONTEXTS):
        ax = axes[row_idx][col_idx]
        mat = mats[ctx].astype(float)
        # Mask diagonal for vmin calculation
        off_diag = mat.values[~np.eye(len(SAMPLES), dtype=bool)]
        valid = off_diag[~np.isnan(off_diag)]
        vmin = max(0, valid.min() - 0.05) if len(valid) > 0 else 0

        annot_mat = mat.map(lambda v: f'{v:.3f}' if not np.isnan(v) else 'N/A')
        sns.heatmap(
            mat,
            ax=ax,
            annot=annot_mat,
            fmt='',
            cmap='Blues',
            vmin=vmin,
            vmax=1.0,
            square=True,
            linewidths=0.5,
            cbar_kws={'shrink': 0.8},
            annot_kws={'size': 11},
        )
        ax.set_title(f'{method} – {ctx}', fontsize=11)
        ax.set_xticklabels(SAMPLES, rotation=30, ha='right')
        ax.set_yticklabels(SAMPLES, rotation=0)

plt.tight_layout()
for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'correlation_heatmaps.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")
plt.close()

# ── 3. Chromosome-level PCA ───────────────────────────────────────────────────
# Use mean methylation per chromosome per sample per context as features.
# This avoids the sparsity issue (different genomes → different site sets).
# Result: 3 sample points in a space defined by chr × context mean methylations.

print("\nBuilding chromosome-level feature matrix for PCA...")
CHROM_ORDER = WHEAT_CHROMS + [f'{i}H' for i in range(1, 8)]

feature_rows = {s: {} for s in SAMPLES}

for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    df = pd.read_csv(fpath, sep='\t',
                     usecols=['Chromosome', 'Mv9', 'Igri', '7Hadd'],
                     low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)
    for col in SAMPLES:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    for chrom in CHROM_ORDER:
        sub = df[df['Chromosome'] == chrom]
        for s in SAMPLES:
            key = f'{ctx}_{chrom}'
            val = sub[s].dropna().mean()
            feature_rows[s][key] = val

feature_df = pd.DataFrame(feature_rows).T   # samples × features
feature_df = feature_df.fillna(feature_df.mean())  # impute missing chr means with global mean

print(f"  Feature matrix shape: {feature_df.shape}  (samples × chr×context means)")

# PCA via SVD
X = feature_df.values.astype(float)
X_c = X - X.mean(axis=0)   # centre features
U, S, Vt = np.linalg.svd(X_c, full_matrices=False)
explained = (S ** 2) / (S ** 2).sum() * 100
pc1 = U[:, 0]
pc2 = U[:, 1] if len(S) > 1 else np.zeros(len(SAMPLES))

fig, ax = plt.subplots(figsize=(6, 5))
for i, sample in enumerate(SAMPLES):
    ax.scatter(pc1[i], pc2[i],
               color=SAMPLE_COLORS[sample], s=250, zorder=5,
               edgecolors='black', linewidths=0.8)
    ax.annotate(sample, xy=(pc1[i], pc2[i]),
                xytext=(10, 6), textcoords='offset points', fontsize=12)

ax.axhline(0, color='grey', lw=0.6, ls='--')
ax.axvline(0, color='grey', lw=0.6, ls='--')
ax.set_xlabel(f'PC1 ({explained[0]:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2 ({explained[1]:.1f}% variance)' if len(explained) > 1 else 'PC2', fontsize=12)
ax.set_title(
    'Chromosome-level PCA of Methylation Profiles\n'
    '(features: mean CpG/CHG/CHH per chromosome)',
    fontsize=11
)
ax.grid(alpha=0.3)
patches = [mpatches.Patch(color=SAMPLE_COLORS[s], label=s) for s in SAMPLES]
ax.legend(handles=patches, fontsize=10)

plt.tight_layout()
for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'pca_chromosome_level.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")
plt.close()

print(f"\nPCA variance: PC1={explained[0]:.1f}%  PC2={explained[1]:.1f}%")
print("\n✓ da03 complete.")
