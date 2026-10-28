# Breast Cancer RNA-seq: Differential Expression Analysis

![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) ![HISAT2](https://img.shields.io/badge/aligner-HISAT2-blue) ![DESeq2](https://img.shields.io/badge/DE-PyDESeq2-green)

## Objective

Identify differentially expressed genes and affected biological pathways between tumour (invasive ductal carcinoma, Luminal B subtype) and adjacent normal tissue from a single breast cancer patient.

## Data

- **NCBI SRA BioProject:** [PRJNA432903](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA432903) (Hong et al., *Data in Brief*, 2018, DOI: [10.1016/j.dib.2018.03.079](https://doi.org/10.1016/j.dib.2018.03.079))
- 10 tumour (C0–C9) + 3 adjacent normal (N1–N3) samples
- Illumina HiSeq 2500, paired-end RNA-seq

## Pipeline

```
SRA (prefetch, fasterq-dump) → HISAT2 (GRCh38) → featureCounts → PyDESeq2 → mygene (Ensembl → symbol) → gseapy / Enrichr (GO BP 2021) → plots
```

| Step | Tool | Script |
|---|---|---|
| Download | SRA Toolkit | `scripts/download_resume.sh` |
| Alignment + counting | HISAT2, Subread featureCounts | `scripts/align_count_clean.sh` |
| Count matrix | pandas | `scripts/build_matrix.py` |
| Differential expression | PyDESeq2 | `scripts/run_deseq2.py` |
| Annotation & plots | mygene, matplotlib, seaborn, scikit-learn | `scripts/analyze_and_plot.py` |
| GO enrichment | gseapy (Enrichr) | `scripts/go_enrichment.py` |

## Results

### Differential expression

**8,971 significant genes** (padj < 0.05 and |log2FC| > 1). Top hits include immunoglobulin genes (*IGHV3-43*, *IGLV3-19*, *IGKV1-27*) and T cell / lymphocyte genes (*CD8A*, *TIGIT*, *LAX1*, *SLA2*, *CXCL13*).

| Volcano plot | PCA |
|---|---|
| ![Volcano plot](figures/volcano_plot.png) | ![PCA](figures/pca_plot.png) |

![Heatmap of top 30 genes](figures/heatmap_top30.png)

### GO enrichment

- **Up-regulated genes:** strong enrichment for immune response, T cell activation and proliferation, interferon-gamma signalling and B cell receptor signalling (adj. p < 1e-5)
- **Down-regulated genes:** only *extracellular matrix assembly* was significant (adj. p = 0.0019); other terms were not significant

### Comparison with the original study

These results independently agree with the *lymphocyte activation* finding reported in the original publication, which supports that the pipeline works correctly.

## Repository Structure

```
rnaseq-meme-kanseri-analizi/
├── raw_data/   # SRR list, run info, metadata (FASTQ files not included)
├── scripts/    # all analysis scripts
├── results/    # count matrix, DESeq2 results, GO tables, alignment summaries
└── figures/    # volcano plot, PCA, heatmap
```

## Limitations

- All samples come from a single patient, so inter-individual variation cannot be controlled and results do not generalise.
- The original paper reports low quality for sample C3 (TIN score 59.6); this analysis kept C3.
- No clear, consistent pathway emerged among down-regulated genes, possibly due to sample size or effect size.
- This is a learning/portfolio project and must not be used for clinical decisions.

## Author

**Atakan Aslan** · Molecular Biology and Genetics, Sivas Cumhuriyet University
