# Analysis pipelines

This folder contains the main **R scripts used during my bachelor's thesis** on transcriptomic and metabolic changes in pancreatic ductal adenocarcinoma (PDAC).

The code is preserved relatively close to the versions used during the research project. It therefore contains exploratory sections, intermediate outputs, extensive comments and some research-specific file paths.

It is intended primarily to document **how the analysis was performed and how the workflow developed**, rather than to function as a polished software package.

For the final methodology, results and interpretation, see the [main repository README](../README.md) and the [full thesis](../Jakob_Urh_Thesis.pdf).

---

## Scripts

### `run_DESeq2.R`

Performs the differential-expression analysis using **DESeq2**.

Main steps include:

- combining RNA-seq count data from CPTAC-3 and TCGA-PAAD;
- matching genes shared between the datasets;
- removing Ensembl version suffixes;
- aggregating duplicated gene identifiers;
- constructing the DESeq2 dataset;
- modelling dataset origin and tumour/normal status;
- filtering low-count genes;
- estimating differential expression;
- shrinking log2 fold changes with `apeglm`;
- identifying genes using adjusted p-value and fold-change thresholds;
- annotating Ensembl identifiers using `biomaRt`;
- exporting intermediate and final DESeq2 results.

The analysis used:

```text
padj < 0.01
|log2FC| >= 1
```

for the final differential-expression candidate set.

---

### `logTPM_analysis.R`

Performs the analysis of **log2-transformed TPM values** for the paired CPTAC-3 cohort.

Main steps include:

- cleaning Ensembl gene identifiers;
- averaging duplicated TPM entries;
- retaining patients with both tumour and normal samples;
- averaging replicate samples where necessary;
- transforming expression using `log2(TPM + 1)`;
- connecting genes to Human1 metabolic subsystems;
- calculating normalized tumour-normal Euclidean distances;
- ranking metabolic subsystems;
- performing paired Wilcoxon tests;
- applying multiple-testing correction;
- calculating ROC-AUC values;
- identifying transcriptomic features with strong tumour-normal discrimination.

This script also contains exploratory visualizations and intermediate analyses used while developing the final methodology.

---

### `flux_analysis_v2.R`

Processes and analyses the **estimated metabolic-flux data**.

The original flux estimates were generated in previous work using context-specific Human1 models and were used as input for this thesis.

Main steps include:

- reading individual sample-level flux files;
- matching flux files to patient metadata;
- constructing a common reaction-by-sample matrix;
- retaining paired tumour and normal patients;
- mapping reactions to Human1 metabolic subsystems;
- calculating normalized tumour-normal Euclidean distances;
- performing paired Wilcoxon tests;
- applying multiple-testing correction;
- calculating reaction-level ROC-AUC values;
- ranking altered reactions and metabolic subsystems.

> **Note:** This script still contains a local Windows path for the original flux-data directory and therefore requires path modification before reuse.

---

### `DESeq2_log2TPM_flux.R`

Integrates the results from the three main analytical layers:

- DESeq2 differential expression;
- TPM-based analysis;
- metabolic-flux analysis.

The script uses Human1 gene-reaction-subsystem relationships to construct a combined table containing information such as:

- gene identifier;
- reaction identifier;
- metabolic subsystem;
- DESeq2 adjusted p-value;
- log2 fold change;
- normalized mean counts;
- TPM AUC;
- TPM paired-Wilcoxon adjusted p-value;
- flux AUC;
- flux paired-Wilcoxon adjusted p-value;
- subsystem rankings.

This integrated table was used to identify candidate features supported across multiple analytical approaches and later served as input for the LASSO analysis.

---

### `ML_v2.R`

Contains the final **LASSO feature-selection and classification workflow**.

The script:

- selects candidate transcriptomic and reaction features;
- uses the same 52 paired patients across analyses;
- splits patients rather than individual samples;
- keeps tumour and normal samples from the same patient in the same partition;
- uses 35 patients for training and 17 for testing;
- performs patient-level internal cross-validation;
- uses `lambda.1se` for feature selection;
- repeats the procedure across 100 random seeds;
- measures feature-selection frequency;
- retains features selected in more than 70% of runs;
- evaluates gene-only, reaction-only and combined models;
- calculates AUC and classification accuracy;
- records model stability and non-convergence.

The final comparison included:

| Model | Features | Median AUC | Median accuracy |
|---|---:|---:|---:|
| Gene | 3 | **0.990** | **0.941** |
| Reaction | 7 | 0.948 | 0.882 |
| Combined | 10 | 0.933 | 0.912 |

The combined model was substantially less stable, with non-convergence in 83% of repeated patient splits.

For interpretation of these values and the feature-selection leakage limitation, see the main thesis README.

---

### `graphs.R`

Contains code used to generate and refine figures for the final thesis and presentation.

The script uses outputs from the analytical pipelines to create publication-style visualizations of results such as:

- distributions of normalized Euclidean distances;
- metabolic subsystem comparisons;
- transcriptomic results;
- flux results;
- integrated candidate features;
- model-performance summaries.

Some labels are written in Slovenian because these figures were created for the final Slovenian-language thesis.

---

## Approximate workflow

The scripts broadly follow this order:

```text
1. run_DESeq2.R
2. logTPM_analysis.R
3. flux_analysis_v2.R
4. DESeq2_log2TPM_flux.R
5. ML_v2.R
6. graphs.R
```

The three initial analyses are partly independent and produce intermediate files that are later combined by `DESeq2_log2TPM_flux.R`.

`ML_v2.R` uses the integrated results together with the cleaned paired transcriptomic and flux datasets.

---

## Important note on reproducibility

This folder is **not currently a self-contained one-command pipeline**.

Several required input files are not redistributed in this repository, including the original transcriptomic data, Human1 model file and metabolic-flux input data.

Some scripts also retain:

- original local file paths;
- intermediate file reads and writes;
- package-installation commands;
- exploratory code;
- plotting experiments;
- comments written while I was learning and developing the analysis.

Reproducing the full workflow therefore requires obtaining the original datasets and adapting file paths to the local environment.

The purpose of preserving these scripts is to provide a transparent record of the computational work behind the thesis rather than to present the project as production software.

---

## Why the original code is preserved

I considered rewriting the scripts into a cleaner pipeline after completing the thesis.

Instead, I chose to preserve the research code relatively close to the versions I actually worked with.

This makes some parts less elegant, but it also documents the development of the project:

- early exploration;
- data-cleaning decisions;
- intermediate checks;
- debugging;
- methodological changes;
- visualization experiments;
- increasingly structured analysis code.

In particular, the later `ML_v2.R` pipeline is considerably more structured than some of the earlier exploratory scripts.

For me, that progression is part of what this repository is intended to show.

---

## Language

Most code comments are written in **Slovenian**, as the scripts were developed alongside a Slovenian-language bachelor's thesis.

Variable and function names therefore contain a mixture of Slovenian and English terminology.
