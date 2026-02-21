#!/usr/bin/env python3
"""
da06_replicate_validation.py
----------------------------
Biological replicate reproducibility analysis for Scientific Data Technical Validation.

Nine samples: Mv9_1/2/3 (wheat), Igri_1/2/3 (barley), 7Hadd_1/2/3 (addition line).

Chromosome filtering per genotype (applied before all analyses):
  - Mv9     : wheat chromosomes only (1A–7D) — Igri has no wheat reads
  - Igri    : barley chromosomes only (1H–7H) — Mv9 has no barley reads
  - 7Hadd   : wheat chromosomes + 7H only  (1H–6H rows are cross-mapping artefacts)

NaN handling:
  - Correlations : pairwise complete cases (dropna on both columns)
  - PCA          : chromosome-level mean methylation per sample (avoids site-level
                   sparsity; missing chromosome means imputed with grand mean)

Outputs (all written to this script's directory):
  replicate_correlation_summary.csv      — Pearson & Spearman r for all 36 sample pairs
  replicate_correlation_heatmap.png/.svg — 9×9 annotated heatmap, genotype colour bands
  replicate_pca.png/.svg                 — chromosome-level PCA, 9 sample points
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from itertools import combinations, product

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = SCRIPT_DIR          # replicate files are in this directory
OUT_DIR    = SCRIPT_DIR

CONTEXTS = ['CpG', 'CHG', 'CHH']

WHEAT_CHROMS  = [f'{i}{s}' for i in range(1, 8) for s in ['A', 'B', 'D']]
BARLEY_CHROMS = [f'{i}H' for i in range(1, 8)]
BARLEY_7H     = ['7H']
CHROM_ORDER   = WHEAT_CHROMS + BARLEY_CHROMS

# Genotype → allowed chromosomes
GENO_CHROMS = {
    'Mv9':    WHEAT_CHROMS,
    'Igri':   BARLEY_CHROMS,
    '7Hadd':  WHEAT_CHROMS + BARLEY_7H,
}

# Sample display order and colours
SAMPLES = ['Mv9_1','Mv9_2','Mv9_3','Igri_1','Igri_2','Igri_3','7Hadd_1','7Hadd_2','7Hadd_3']
GENO_COLOR = {
    'Mv9':   '#2166AC',
    'Igri':  '#D6604D',
    '7Hadd': '#4DAC26',
}
SAMPLE_COLOR = {s: GENO_COLOR[s.rsplit('_',1)[0]] for s in SAMPLES}

def genotype(sample):
    return sample.rsplit('_', 1)[0]

def allowed_chroms(sample):
    return GENO_CHROMS[genotype(sample)]

def shared_chroms(s1, s2):
    """Return chromosomes present in both samples' allowed sets."""
    return list(set(allowed_chroms(s1)) & set(allowed_chroms(s2)))

def short_name(col):
    """
    Handles both CpG format:  '7Hadd_1_S7_R1_...'
    and CHG/CHH format:       'CHG_context_7Hadd_1_S7_R1_...'
    → always returns e.g. '7Hadd_1'
    """
    import re
    m = re.search(r'(7Hadd|Igri|Mv9)_([123])', col)
    if m:
        return f'{m.group(1)}_{m.group(2)}'
    return col

print("=== da06: Replicate Reproducibility Validation ===")

# ── 1. Load, rename, filter all contexts ────────────────────────────────────

data = {}   # data[ctx] = DataFrame with short sample names + Chromosome

for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'{ctx}_all_replicates.txt')
    print(f"\nLoading {ctx}...")
    df = pd.read_csv(fpath, sep='\t', low_memory=False)
    df['Chromosome'] = df['Chromosome'].astype(str)

    # Rename long sample columns to short names
    rename = {c: short_name(c)
              for c in df.columns
              if any(g in c for g in ['7Hadd','Igri','Mv9'])}
    df = df.rename(columns=rename)

    for s in SAMPLES:
        df[s] = pd.to_numeric(df[s], errors='coerce')

    data[ctx] = df[['Chromosome'] + SAMPLES]
    print(f"  {len(df):,} probes, columns: {SAMPLES}")

# ── 2. Pairwise correlations ─────────────────────────────────────────────────

print("\nComputing pairwise correlations...")
corr_rows = []

# Build 9×9 matrices for heatmap (average Pearson across contexts)
pearson_mat_ctx  = {}   # per-context
avg_pearson      = pd.DataFrame(np.nan, index=SAMPLES, columns=SAMPLES)

for ctx in CONTEXTS:
    df = data[ctx]
    pm = pd.DataFrame(np.nan, index=SAMPLES, columns=SAMPLES)

    for s in SAMPLES:
        pm.loc[s, s] = 1.0

    for s1, s2 in combinations(SAMPLES, 2):
        chroms = shared_chroms(s1, s2)
        sub = df[df['Chromosome'].isin(chroms)][[s1, s2]].dropna()
        n = len(sub)
        if n < 20:
            print(f"  SKIP {s1}↔{s2} ({ctx}): only {n} shared sites")
            continue

        r_p, _ = stats.pearsonr(sub[s1], sub[s2])
        r_s, _ = stats.spearmanr(sub[s1], sub[s2])
        pm.loc[s1, s2] = pm.loc[s2, s1] = round(r_p, 4)

        within = genotype(s1) == genotype(s2)
        corr_rows.append({
            'Context':        ctx,
            'Sample_A':       s1,
            'Sample_B':       s2,
            'Genotype_A':     genotype(s1),
            'Genotype_B':     genotype(s2),
            'Comparison_type':'within' if within else 'between',
            'Chromosomes':    'wheat+7H' if set(chroms) == set(WHEAT_CHROMS+BARLEY_7H)
                              else ('barley' if chroms == BARLEY_CHROMS
                              else ('7H' if chroms == BARLEY_7H else 'wheat')),
            'N_shared_sites': n,
            'Pearson_r':      round(r_p, 4),
            'Spearman_r':     round(r_s, 4),
        })

    pearson_mat_ctx[ctx] = pm

# Average Pearson across contexts for summary heatmap
stack = np.stack([pearson_mat_ctx[c].values.astype(float) for c in CONTEXTS], axis=0)
avg_pearson = pd.DataFrame(
    np.nanmean(stack, axis=0),
    index=SAMPLES, columns=SAMPLES
)

corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(os.path.join(OUT_DIR, 'replicate_correlation_summary.csv'), index=False)

# Print within-genotype summary
print("\n── Within-genotype Pearson r (mean ± std across contexts) ──")
for geno in ['Mv9','Igri','7Hadd']:
    sub = corr_df[(corr_df['Comparison_type']=='within') &
                  (corr_df['Genotype_A']==geno)]
    print(f"  {geno}: mean r={sub['Pearson_r'].mean():.4f}  "
          f"std={sub['Pearson_r'].std():.4f}  "
          f"min={sub['Pearson_r'].min():.4f}")

# ── 3. Correlation heatmap (9×9, all contexts) ───────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(20, 6))
fig.suptitle('Pairwise Pearson Correlations Across Biological Replicates',
             fontsize=14, fontweight='bold')

for ax, ctx in zip(axes, CONTEXTS):
    mat = pearson_mat_ctx[ctx].astype(float)
    off = mat.values[~np.eye(9, dtype=bool)]
    vmin = max(0, np.nanmin(off) - 0.02)

    annot = mat.map(lambda v: f'{v:.3f}' if not np.isnan(v) else '')
    g = sns.heatmap(mat, ax=ax, annot=annot, fmt='', cmap='Blues',
                    vmin=vmin, vmax=1.0, square=True,
                    linewidths=0.4, linecolor='white',
                    cbar_kws={'shrink': 0.7},
                    annot_kws={'size': 7})

    # Colour-code tick labels by genotype
    ax.set_xticklabels(SAMPLES, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(SAMPLES, rotation=0, fontsize=8)
    for tick, sample in zip(ax.get_xticklabels(), SAMPLES):
        tick.set_color(GENO_COLOR[genotype(sample)])
    for tick, sample in zip(ax.get_yticklabels(), SAMPLES):
        tick.set_color(GENO_COLOR[genotype(sample)])

    # Draw boxes around within-genotype blocks (3×3 each)
    for start in [0, 3, 6]:
        ax.add_patch(plt.Rectangle((start, start), 3, 3,
                     fill=False, edgecolor='black', lw=2.0, zorder=5))

    ax.set_title(f'{ctx}', fontsize=12)

# Legend
patches = [mpatches.Patch(color=c, label=g) for g, c in GENO_COLOR.items()]
fig.legend(handles=patches, loc='lower center', ncol=3,
           fontsize=10, bbox_to_anchor=(0.5, -0.04))

plt.tight_layout(rect=[0, 0.04, 1, 1])
for ext in ['png', 'svg']:
    path = os.path.join(OUT_DIR, f'replicate_correlation_heatmap.{ext}')
    plt.savefig(path, dpi=180, bbox_inches='tight')
    print(f"Saved: {path}")
plt.close()

# ── 4. Chromosome-level PCA (9 samples) ──────────────────────────────────────

print("\nBuilding chromosome-level feature matrix for PCA (9 samples)...")
feature_rows = {s: {} for s in SAMPLES}

for ctx in CONTEXTS:
    df = data[ctx]
    for chrom in CHROM_ORDER:
        sub = df[df['Chromosome'] == chrom]
        for s in SAMPLES:
            key = f'{ctx}_{chrom}'
            # Only use chromosome if it's in the sample's allowed set
            if chrom in allowed_chroms(s):
                val = sub[s].dropna().mean()
            else:
                val = np.nan
            feature_rows[s][key] = val

feat_df = pd.DataFrame(feature_rows).T   # 9 samples × 84 features
# Impute missing feature means column-wise (missing = chromosome not in genome)
feat_df = feat_df.apply(lambda col: col.fillna(col.mean()), axis=0)
# Drop features with all-NaN (shouldn't occur after imputation)
feat_df = feat_df.dropna(axis=1)

print(f"  Feature matrix: {feat_df.shape}  (samples × chr×context means)")

X = feat_df.values.astype(float)
X_c = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(X_c, full_matrices=False)
explained = (S**2) / (S**2).sum() * 100
pc1, pc2 = U[:, 0], U[:, 1]

fig, ax = plt.subplots(figsize=(7, 6))

for i, sample in enumerate(SAMPLES):
    geno = genotype(sample)
    ax.scatter(pc1[i], pc2[i],
               color=GENO_COLOR[geno], s=220, zorder=5,
               edgecolors='black', linewidths=0.8)
    ax.annotate(sample, xy=(pc1[i], pc2[i]),
                xytext=(8, 5), textcoords='offset points', fontsize=9)

# Draw convex hulls / ellipses around each genotype's 3 points
from matplotlib.patches import Ellipse
for geno, color in GENO_COLOR.items():
    idx = [i for i, s in enumerate(SAMPLES) if genotype(s) == geno]
    pts = np.column_stack([pc1[idx], pc2[idx]])
    if len(pts) >= 3:
        center = pts.mean(axis=0)
        # Simple bounding circle
        r = np.max(np.linalg.norm(pts - center, axis=1)) * 1.4
        circle = plt.Circle(center, r, color=color, alpha=0.10, zorder=1)
        ax.add_patch(circle)

ax.axhline(0, color='grey', lw=0.6, ls='--')
ax.axvline(0, color='grey', lw=0.6, ls='--')
ax.set_xlabel(f'PC1 ({explained[0]:.1f}% variance)', fontsize=12)
ax.set_ylabel(f'PC2 ({explained[1]:.1f}% variance)', fontsize=12)
ax.set_title(
    'PCA of Biological Replicates (chromosome-level features)\n'
    'Shaded areas: within-genotype clusters',
    fontsize=11
)
ax.grid(alpha=0.3)
patches = [mpatches.Patch(color=c, label=g) for g, c in GENO_COLOR.items()]
ax.legend(handles=patches, fontsize=10)

plt.tight_layout()
for ext in ['png', 'svg']:
    path = os.path.join(OUT_DIR, f'replicate_pca.{ext}')
    plt.savefig(path, dpi=180, bbox_inches='tight')
    print(f"Saved: {path}")
plt.close()

print(f"\nPCA variance: PC1={explained[0]:.1f}%  PC2={explained[1]:.1f}%  PC3={explained[2]:.1f}%")

# ── 5. Print summary for manuscript ──────────────────────────────────────────

print("\n══ Summary for Technical Validation ══")
for ctx in CONTEXTS:
    print(f"\n{ctx}:")
    for comp_type in ['within', 'between']:
        sub = corr_df[(corr_df['Context']==ctx) &
                      (corr_df['Comparison_type']==comp_type)]
        if sub.empty: continue
        label = 'Within-genotype ' if comp_type=='within' else 'Between-genotype'
        print(f"  {label}: mean r={sub['Pearson_r'].mean():.4f}  "
              f"range [{sub['Pearson_r'].min():.4f}–{sub['Pearson_r'].max():.4f}]  "
              f"n={len(sub)} pairs")

print("\n✓ da06 complete. Outputs written to:", OUT_DIR)
