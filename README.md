[![DOI](https://zenodo.org/badge/1162721573.svg)](https://doi.org/10.5281/zenodo.18715140)

# RRBS Data Quality & Technical Validation Scripts

**Associated manuscript:** *"Reduced Representation Bisulfite Sequencing (RRBS) dataset of a wheat–barley 7H chromosome addition line"* — submitted to *Scientific Data*.

These Python scripts were used to produce all Technical Validation figures and supporting tables in the manuscript. They assess data quality, reproducibility, and genomic coverage of an RRBS dataset generated from three plant genotypes: *Triticum aestivum* cv. Mv9 (wheat control), *Hordeum vulgare* cv. Igri (barley control), and the wheat–barley 7H chromosome addition line (7Hadd).

---

## Repository structure

```
data_analysis/
├── da01_data_quality.py        # Site counts, methylation statistics, bisulfite conversion
├── da02_distributions.py       # Methylation distribution plots (KDE) per context
├── da03_correlation_pca.py     # Pairwise sample correlations and chromosome-level PCA
├── da04_chromosome_coverage.py # Probe coverage per chromosome per sample
├── da05_genomic_features.py    # Genomic feature distribution of differentially methylated probes
└── README.md

# Upstream preprocessing scripts (required to reproduce input files from raw data):
../step6_filter_context_files.py   # Filters raw SeqMonk exports (≥2 non-missing samples)
../step8_smart_clean.py            # Removes all-NaN rows, preserves biological zeros → clean_*.txt
../step4_add_status_columns.py     # Adds Hyper/Hypo status columns to DMP table
```

---

## Input files

The five `da*.py` scripts in this directory require the following pre-processed input files, which must be placed in the **parent directory** (`../`) relative to this folder.

### `clean_CpG_all.txt`, `clean_CHG_all.txt`, `clean_CHH_all.txt`
*Used by: da01, da02, da03, da04*

Tab-separated files containing genome-wide RRBS methylation values (%) for all three samples. These are the **primary deposited datasets**, available from the data repository (see Data Records in the manuscript).

**Column structure:**
| Column | Description |
|--------|-------------|
| `Chromosome` | Chromosome identifier (e.g. `1A`, `7H`) |
| `Start` | Probe start position (0-based) |
| `End` | Probe end position |
| `Mv9` | Methylation % in *T. aestivum* cv. Mv9 |
| `Igri` | Methylation % in *H. vulgare* cv. Igri |
| `7Hadd` | Methylation % in the 7H addition line |

**Provenance:** Raw SeqMonk exports (`CpG_all.txt`, `CHG_all.txt`, `CHH_all.txt`) were processed by:
1. `step6_filter_context_files.py` — keeps rows with data in ≥2 samples
2. `step8_smart_clean.py` — converts non-numeric values to NaN, drops rows where **all three** samples are NaN (preserving biological zeros)

### `filtered_significant_methylation.csv`
*Used by: da05*

Comma-separated file of statistically significant differentially methylated probes (DMPs) identified by SeqMonk logistic regression (p < 0.05 after Benjamini–Hochberg correction, minimum 10 observations). Contains columns for chromosome, position, methylation context, genomic region annotation, comparison group (`Wheat_Genom` or `Barley_7H`), and methylation direction (`Mv9_7Hadd_Status`, `Igri_7Hadd_Status`: `hyper`/`hypo`).

**Provenance:** SeqMonk statistical output annotated with genomic region and Hyper/Hypo status by `step4_add_status_columns.py`.

---

## Requirements

**Python ≥ 3.10** with the following packages:

| Package | Version tested | Purpose |
|---------|---------------|---------|
| `pandas` | ≥2.0 | Data loading and manipulation |
| `numpy` | ≥1.24 | Numerical operations, SVD for PCA |
| `matplotlib` | ≥3.7 | Plotting backend |
| `seaborn` | ≥0.12 | KDE plots, heatmaps |
| `scipy` | ≥1.10 | Pearson and Spearman correlations |

> **Note:** `pandas ≥ 2.0` is required; `DataFrame.applymap()` has been removed — the scripts use `DataFrame.map()` instead.

### Installation via conda (recommended)

A `conda` environment file is provided in the repository root:

```bash
conda env create -f ../environment.yaml
conda activate gemini_rrbs
```

Or install manually:

```bash
pip install pandas>=2.0 numpy matplotlib seaborn scipy
```

---

## Running the scripts

All scripts must be run **from the `data_analysis/` directory**. They resolve input paths relative to their own location (`../clean_CpG_all.txt`, etc.) and write all outputs to the same directory.

```bash
cd data_analysis/

# Option A: with conda environment active
python da01_data_quality.py
python da02_distributions.py
python da03_correlation_pca.py
python da04_chromosome_coverage.py
python da05_genomic_features.py

# Option B: direct path to interpreter (no activation needed)
/path/to/envs/gemini_rrbs/bin/python da01_data_quality.py
```

Scripts are independent and can be run in any order. No script modifies the input files.

---

## Script descriptions and outputs

### `da01_data_quality.py`
Calculates per-sample site counts, methylation summary statistics (mean, median, SD, range), data completeness (sites present in 1/2/3 samples), and estimates bisulfite conversion efficiency from CHH methylation on wheat chromosomes.

| Output file | Description |
|-------------|-------------|
| `data_quality_summary.csv` | Site counts and methylation stats per context × sample |
| `data_completeness.csv` | Sites shared across 1, 2, or all 3 samples per context |
| `bisulfite_conversion_estimate.csv` | Mean CHH methylation and conversion efficiency per sample |
| `per_chromosome_site_counts.csv` | Site counts per chromosome × context × sample |

---

### `da02_distributions.py`
Generates a three-panel KDE figure showing the full methylation value distribution for each context (CpG, CHG, CHH) and all three samples. CpG distributions display the expected bimodal pattern; CHH is right-skewed and low, confirming data quality.

| Output file | Description |
|-------------|-------------|
| `methylation_distributions.png` | Figure 2 (manuscript) — raster format, 180 dpi |
| `methylation_distributions.svg` | Vector format for publication |

---

### `da03_correlation_pca.py`
Computes pairwise Pearson and Spearman correlations using biologically appropriate chromosome subsets (Mv9 ↔ 7Hadd on wheat chromosomes; Igri ↔ 7Hadd on chromosome 7H), and plots correlation heatmaps. Also performs chromosome-level PCA using mean methylation per chromosome per sample per context as features (avoids the sparse-overlap problem inherent to comparing different genomes).

| Output file | Description |
|-------------|-------------|
| `correlation_summary.csv` | Pearson r, Spearman r, p-values, N per context × pair |
| `correlation_heatmaps.png` | Figure 4 (manuscript) — 3×2 heatmap grid |
| `correlation_heatmaps.svg` | Vector format |
| `pca_chromosome_level.png` | Figure 5 (manuscript) — chromosome-level PCA |
| `pca_chromosome_level.svg` | Vector format |

---

### `da04_chromosome_coverage.py`
Plots the number of methylation sites with data per chromosome per sample for each context as a grouped bar chart. Separates wheat and barley chromosomes with a vertical divider.

| Output file | Description |
|-------------|-------------|
| `chromosome_coverage.png` | Figure 3 (manuscript) — raster format |
| `chromosome_coverage.svg` | Vector format |
| `chromosome_coverage_detail.csv` | Site counts per chromosome × context × sample |

---

### `da05_genomic_features.py`
Analyses the distribution of significant DMPs across annotated genomic regions (promoter / gene body / downstream) for wheat and barley 7H comparisons, and breaks down hyper- vs. hypomethylation per region.

| Output file | Description |
|-------------|-------------|
| `genomic_features_summary.csv` | DMP counts and proportions per genome × context × region |
| `genomic_features_barplot.png` | Stacked bar chart of DMP distribution |
| `genomic_features_barplot.svg` | Vector format |
| `hyper_hypo_breakdown.csv` | Hyper/hypo counts per genome × context × region |
| `hyper_hypo_barplot.png` | Stacked bar chart of hyper/hypo counts |
| `hyper_hypo_barplot.svg` | Vector format |

---

## Sample labels

| Label | Full name | Genome composition |
|-------|-----------|--------------------|
| `Mv9` | *Triticum aestivum* cv. Mv9 | Wheat (AABBDD) |
| `Igri` | *Hordeum vulgare* cv. Igri | Barley (7× 1H–7H) |
| `7Hadd` | Wheat–barley 7H addition line | Wheat + barley chromosome 7H |

**Chromosome naming:** Wheat chromosomes: `1A`, `1B`, `1D` … `7A`, `7B`, `7D`. Barley chromosomes: `1H`–`7H`. In the 7H addition line, only `7H` is a genuine alien chromosome; `1H`–`6H` rows originate from cross-mapping artefacts in the Igri barley control.

---

## License

MIT License. See `LICENSE` file.

---

## Citation

If you use these scripts, please cite:

> [Author list] (*year*). Reduced Representation Bisulfite Sequencing (RRBS) dataset of a wheat–barley 7H chromosome addition line. *Scientific Data*. DOI: [to be assigned]
