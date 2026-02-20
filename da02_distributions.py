#!/usr/bin/env python3
"""
da02_distributions.py
---------------------
Methylation distribution visualisation for Scientific Data submission.

Outputs:
  methylation_distributions.png / .svg
      3-panel figure: one row per context (CpG / CHG / CHH).
      Each panel shows KDE + rug for all three samples.
      CpG bimodality is expected; CHH should be right-skewed and low.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..')
OUT_DIR    = SCRIPT_DIR

CONTEXTS = ['CpG', 'CHG', 'CHH']
SAMPLES  = ['Mv9', 'Igri', '7Hadd']

SAMPLE_COLORS = {
    'Mv9':    '#2166AC',   # blue  – wheat control
    'Igri':   '#D6604D',   # red   – barley control
    '7Hadd':  '#4DAC26',   # green – addition line
}

SAMPLE_LABELS = {
    'Mv9':    'Mv9 (wheat)',
    'Igri':   'Igri (barley)',
    '7Hadd':  '7H addition line',
}

print("=== da02: Methylation Distributions ===")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('DNA Methylation Distributions by Context', fontsize=15, fontweight='bold', y=1.01)

for ax, ctx in zip(axes, CONTEXTS):
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    df = pd.read_csv(fpath, sep='\t', usecols=['Mv9', 'Igri', '7Hadd'], low_memory=False)

    for sample in SAMPLES:
        vals = pd.to_numeric(df[sample], errors='coerce').dropna()
        # Downsample for speed (keep max 200k points)
        if len(vals) > 200_000:
            vals = vals.sample(200_000, random_state=42)
        sns.kdeplot(
            vals,
            ax=ax,
            color=SAMPLE_COLORS[sample],
            label=SAMPLE_LABELS[sample],
            linewidth=2.0,
            fill=True,
            alpha=0.18,
        )

    ax.set_title(f'{ctx} methylation', fontsize=13)
    ax.set_xlabel('Methylation (%)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_xlim(-2, 102)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # For CHH: zoom x-axis to 0-30 % because values are low
    if ctx == 'CHH':
        ax.set_xlim(-0.5, 30)
        ax.set_title('CHH methylation (x-axis: 0–30 %)', fontsize=13)

plt.tight_layout()

for ext in ['png', 'svg']:
    outpath = os.path.join(OUT_DIR, f'methylation_distributions.{ext}')
    plt.savefig(outpath, dpi=180, bbox_inches='tight')
    print(f"Saved: {outpath}")

plt.close()

# ── Per-context, per-sample summary statistics printed for reference ──────────
print("\nSummary statistics:")
for ctx in CONTEXTS:
    fpath = os.path.join(DATA_DIR, f'clean_{ctx}_all.txt')
    df = pd.read_csv(fpath, sep='\t', usecols=['Mv9', 'Igri', '7Hadd'], low_memory=False)
    print(f"\n{ctx}:")
    for sample in SAMPLES:
        v = pd.to_numeric(df[sample], errors='coerce').dropna()
        print(f"  {sample:8s}  n={len(v):>8,}  mean={v.mean():.2f}%  median={v.median():.2f}%  std={v.std():.2f}")

print("\n✓ da02 complete.")
