# Integrative Analysis of Transcriptomic and Metabolic Features in Pancreatic Ductal Adenocarcinoma

Bachelor's thesis project investigating transcriptomic and metabolic changes in **pancreatic ductal adenocarcinoma (PDAC)** and their potential for distinguishing tumour from normal samples.

**BSc Biochemistry — University of Ljubljana, Faculty of Chemistry and Chemical Technology, 2026**

**Mentor:** Prof. Dr. Miha Moškon  
**Working mentor:** Assistant Ana Halužan Vasle

Original thesis title:

> *Integrativna analiza transkriptomskih in metabolnih značilnosti pri pankreatičnem duktalnem adenokarcinomu*

---

## Quick links

- [Full bachelor's thesis](Jakob_Urh_Thesis.pdf)
- [Final thesis defense presentation](Thesis_presentation.pdf)
- [Analysis pipelines and research code](analyses_pipelines/)
- [Progress presentations](progress_presentation/)

The thesis and presentations are primarily written in Slovenian. This README provides an English overview of the project, methodology, results and limitations.

---

# Overview

Pancreatic ductal adenocarcinoma (PDAC) is characterized by major changes in cellular metabolism and remains difficult to diagnose and treat.

The aim of this thesis was to investigate whether transcriptomic data and estimated metabolic fluxes could be integrated to identify features associated with PDAC.

Rather than looking only at differential gene expression, I compared several layers of information:

- RNA-seq read counts;
- normalized gene expression using TPM;
- estimated metabolic reaction fluxes;
- metabolic subsystems defined by the Human1 genome-scale metabolic model;
- combinations of transcriptomic and metabolic features selected using LASSO logistic regression.

The project therefore combined **bioinformatics, statistics, systems biology and machine learning**.

---

# Research questions

The project focused on three main questions.

### 1. Do transcriptomic and metabolic analyses identify the same altered metabolic subsystems?

Metabolic subsystems were ranked according to normalized Euclidean distances between paired tumour and normal samples using:

- log2-transformed TPM values;
- estimated metabolic reaction fluxes.

### 2. Can multiple analytical approaches be integrated to identify stronger candidate features?

Results from:

- DESeq2;
- TPM analysis;
- paired Wilcoxon testing;
- ROC-AUC analysis;
- metabolic flux analysis;

were combined through gene-reaction relationships defined in Human1.

### 3. Can combinations of features classify tumour and normal samples?

LASSO logistic regression was used to construct:

- a gene-based model;
- a reaction-based model;
- a combined gene + reaction model.

---

# Hypotheses

Three hypotheses were defined before the analysis.

| Hypothesis | Outcome |
|---|---|
| The most altered metabolic subsystems identified using transcriptomic and flux data would largely overlap. | **Not supported** |
| Integrating DESeq2, TPM and metabolic-flux results would identify promising features for tumour-normal separation. | **Supported** |
| A combined gene + reaction LASSO model would be more predictive and stable than separate models. | **Not supported** |

One of the useful outcomes of the project was therefore that **two of the three original hypotheses were rejected**.

The analyses did not simply confirm the initial expectations.

---

# Data

## CPTAC-3

The main dataset contained bulk RNA-seq data from the CPTAC-3 pancreatic cancer cohort.

The original dataset included:

- **179 tumour samples**
- **66 normal samples**
- **168 patients**

For paired analyses, only patients containing at least one tumour and one normal sample were retained.

Where multiple technical or biological replicates existed for the same patient and sample type, their values were averaged.

The final paired dataset contained:

- **52 patients**
- **104 samples**
- one tumour and one normal observation per patient after aggregation.

This paired design was used for:

- TPM analysis;
- metabolic-flux analysis;
- subsystem comparison;
- paired Wilcoxon testing;
- LASSO modelling.

---

## TCGA-PAAD

TCGA-PAAD data were additionally included in the **DESeq2 differential-expression analysis**.

The dataset contained:

- **179 tumour samples**
- **4 normal samples**

Including both CPTAC-3 and TCGA-PAAD increased the number of samples available for estimation of differential expression.

TCGA-PAAD was not used for the paired TPM, flux or LASSO analyses.

---

# Analysis workflow

A simplified representation of the project is:

```text
                         RNA-seq data
                              │
                 ┌────────────┴────────────┐
                 │                         │
            Raw counts                    TPM
                 │                         │
              DESeq2                  log2(TPM + 1)
                 │                         │
                 │                  Human1 annotation
                 │                         │
                 │                 Subsystem analysis
                 │                         │
                 │                 Wilcoxon + AUC
                 │                         │
                 │                         │
                 │              Estimated metabolic fluxes
                 │                         │
                 │                  Human1 subsystems
                 │                         │
                 │                 Wilcoxon + AUC
                 │                         │
                 └────────────┬────────────┘
                              │
                     Integration of results
                              │
                   Gene-reaction candidates
                              │
                    LASSO feature selection
                              │
           ┌──────────────────┼──────────────────┐
           │                  │                  │
       Gene model       Reaction model     Combined model
```

---

# Data preprocessing

The original CPTAC data contained approximately 60,660 genes identified using Ensembl identifiers.

Version suffixes were removed from Ensembl IDs.

This produced a small number of duplicate identifiers.

These duplicates were treated differently depending on the data type:

- RNA-seq counts were **summed**, because read counts are additive;
- TPM values were **averaged**, because they were already normalized.

For the paired CPTAC analysis, patients without both tumour and normal samples were removed.

After preprocessing, the paired dataset contained:

- **52 patients**
- **104 samples**
- approximately **60,616 genes**

TPM values were transformed using:

```text
log2(TPM + 1)
```

before downstream distance-based analyses.

---

# Human1 metabolic model

The analysis used **Human1 v1.19.0**, a genome-scale reconstruction of human metabolism.

The version used in this project contained approximately:

- **148 metabolic subsystems**
- **2,887 genes**
- **12,971 reactions**

Gene-protein-reaction relationships were used to connect transcriptomic measurements to metabolic reactions and subsystems.

Because individual genes can be connected to multiple reactions or subsystems, the integration process required resolving many-to-many relationships between biological entities.

This became an important part of the analysis.

---

# Transcriptomic analysis

## DESeq2

Differential gene expression was analysed using raw RNA-seq counts with **DESeq2**.

The analysis incorporated data from both CPTAC-3 and TCGA-PAAD.

Genes were considered candidates when:

```text
adjusted p-value < 0.01
|log2 fold change| ≥ 1
```

After filtering:

**6,897 genes** showed statistically significant differential expression.

Of these:

- **2,417** had higher expression in tumour samples;
- **4,480** had lower expression in tumour samples.

The purpose of DESeq2 was primarily to identify statistically supported differences in gene expression.

Statistical significance alone, however, does not necessarily imply strong ability to discriminate between tumour and normal samples.

For this reason, the project also included a separate TPM-based analysis using AUC.

---

# TPM analysis

The TPM analysis focused on metabolic genes associated with Human1.

For each gene:

- tumour and normal expression values were compared using the **paired Wilcoxon test**;
- p-values were corrected using the **Benjamini-Hochberg procedure**;
- an **ROC-AUC** value was calculated to measure tumour-normal discrimination.

Because the direction of discrimination was not the main interest, AUC values below 0.5 were transformed using:

```text
1 - AUC
```

resulting in AUC values between 0.5 and 1.

Among approximately **2,884 metabolic genes**:

- **1,523** had adjusted Wilcoxon p-values < 0.01;
- **402** had AUC > 0.80;
- **58** had AUC > 0.90.

The highest observed gene-level AUC was:

**SLC2A1 — AUC = 0.976**

No individual gene achieved perfect separation between tumour and normal samples.

---

# Metabolic subsystem analysis

I wanted to compare broader metabolic changes rather than only individual genes or reactions.

For every patient and metabolic subsystem, a normalized Euclidean distance was calculated between the tumour and normal sample.

For transcriptomic data this was based on:

```text
log2(TPM + 1)
```

while the flux analysis used estimated reaction fluxes.

The distance was normalized according to the number of genes or reactions in the subsystem.

Because subsystem-distance distributions were often asymmetric, **median distance** was used for ranking rather than the mean.

The interquartile range was also calculated to examine variability between patients.

---

## Transcriptomic subsystem changes

Among the most altered subsystems according to TPM were processes related to:

- cholesterol ester synthesis and hydrolysis;
- acylglycerol metabolism;
- retinol metabolism;
- protein assembly;
- protein modification;
- vitamin A metabolism;
- starch and sucrose metabolism;
- protein degradation.

---

## Flux-based subsystem changes

The most altered subsystems according to estimated metabolic fluxes included:

- pentose phosphate pathway;
- glycolysis / gluconeogenesis;
- pyrimidine metabolism;
- pentose and glucuronate interconversions;
- arginine and proline metabolism;
- alanine, aspartate and glutamate metabolism;
- glycine, serine and threonine metabolism;
- citric-acid-cycle-related reactions.

The highest-ranked transcriptomic and flux-based subsystems showed very little overlap.

This caused the **first hypothesis to be rejected**.

It also suggested that transcript abundance and model-estimated fluxes capture different layers of cellular metabolism.

---

# Metabolic flux analysis

Estimated metabolic fluxes were analysed using a workflow analogous to the TPM analysis.

For each reaction:

- paired tumour and normal values were compared;
- paired Wilcoxon tests were performed;
- p-values were corrected using Benjamini-Hochberg;
- AUC was calculated.

Approximately **9,760 unique reactions** were analysed.

Of these:

- **462 reactions** had adjusted Wilcoxon p-values < 0.01;
- only **8 reactions** had AUC > 0.80.

This was substantially fewer highly discriminative features than observed at the transcriptomic level.

The results also demonstrated an important distinction between:

**statistical significance**

and

**practical discriminatory ability**.

Many reactions differed statistically between tumour and normal samples, but relatively few strongly separated the two groups.

---

# Origin of the metabolic-flux data

The metabolic-flux estimates used in this thesis were obtained from previous work by **Filip Petrovič**.

In that work, transcriptomic data were used to construct context-specific genome-scale metabolic models using:

- Human1;
- ftINIT;
- flux balance analysis (FBA).

My thesis did **not** reconstruct those context-specific models from scratch.

Instead, I used the resulting flux estimates for:

- paired tumour-normal comparison;
- subsystem analysis;
- statistical testing;
- AUC calculation;
- integration with transcriptomic data;
- feature selection and LASSO modelling.

This distinction is important for understanding the scope of the project.

---

# Integration of transcriptomic and metabolic results

Results from the three major analyses were joined using:

- Ensembl gene identifiers;
- Human1 reaction identifiers;
- gene-protein-reaction relationships;
- metabolic subsystem annotations.

The integrated table combined information including:

- DESeq2 adjusted p-values;
- log2 fold changes;
- TPM AUC;
- TPM Wilcoxon adjusted p-values;
- flux AUC;
- flux Wilcoxon adjusted p-values;
- Human1 reactions;
- metabolic subsystems.

As increasingly strict criteria were combined, the candidate space became much smaller.

| Analysis | Candidate features |
|---|---:|
| DESeq2 | 6,897 genes |
| TPM, AUC ≥ 0.80 | 402 genes |
| Flux, AUC ≥ 0.80 | 8 reactions |
| DESeq2 + TPM | 94 genes |
| DESeq2 + TPM + flux | **5 genes / 4 reactions** |

---

# Integrated candidate features

The final integration identified **five genes connected to four metabolic reactions**.

| Gene | Human1 reaction | Interpretation | LFC | TPM AUC | Flux AUC |
|---|---|---|---:|---:|---:|
| **RBP2** | MAR06628 | Retinol transport | -5.75 | 0.803 | 0.806 |
| **STRA6** | MAR06628 | Retinol transport | 1.35 | 0.865 | 0.806 |
| **SLC16A3** | MAR06022 | Acetate / α-ketoglutarate transport | 1.62 | 0.952 | 0.845 |
| **SRD5A3** | MAR01976 | Testosterone → DHT | 1.23 | 0.950 | 0.825 |
| **GPX2** | MAR04116 | H₂O₂ reduction / glutathione oxidation | 1.59 | 0.846 | 0.805 |

These candidates were related to:

- retinol transport;
- monocarboxylate transport;
- androgen metabolism;
- oxidative-stress regulation.

The second hypothesis was therefore **supported**: integrating multiple layers of analysis produced a small set of candidate transcriptomic-metabolic features.

These should be interpreted as **exploratory candidates**, not experimentally validated biomarkers.

---

# Transcript expression and metabolic flux do not always agree

One of the more interesting results was that increased expression of a gene did not necessarily correspond to increased flux through its associated reaction.

For the four integrated reactions, the direction of flux change was also not identical across all 52 patients.

This highlights several biological and modelling issues:

- PDAC is highly heterogeneous;
- bulk RNA-seq measures mixed cell populations;
- transcript abundance does not directly measure enzyme activity;
- metabolic fluxes are model estimates rather than direct experimental measurements;
- several genes can contribute to the same reaction;
- model structure influences gene-reaction relationships.

This was one of the reasons I found the integration step more informative than analysing each layer independently.

---

# Human1 version dependence

The project used **Human1 v1.19.0**.

An important example of model-version dependence involves:

**SLC16A3 — MAR06022**

In Human1 v1.19.0, SLC16A3 was associated with reaction MAR06022.

That association was removed in Human1 v2.0.0.

This means that biological interpretation can depend not only on the experimental data but also on the specific version of the metabolic reconstruction being used.

Re-running the analysis with newer versions of Human1 is therefore an important future improvement.

---

# LASSO logistic regression

The final stage used **L1-regularized logistic regression (LASSO)** to identify combinations of features capable of classifying tumour and normal samples.

Three models were tested:

1. gene features only;
2. metabolic reactions only;
3. combined genes + reactions.

---

## Train/test strategy

The paired dataset contained 52 patients.

For each run:

- **35 patients** were assigned to the training set;
- **17 patients** were assigned to the test set.

Splitting was performed at the **patient level**, ensuring that the tumour and normal samples from the same patient always remained in the same partition.

The process was repeated across **100 random train/test splits**.

Within training data, **10-fold cross-validation** was used to determine the regularization parameter.

The more conservative:

```text
lambda.1se
```

was selected to favour simpler models.

---

# Stability selection

Features were ranked according to how often LASSO selected them across the 100 random splits.

A feature was considered stable when it appeared in at least:

**70 / 100 runs**

This resulted in:

- **3 stable transcriptomic features**
- **7 stable reaction features**

These features were then used in the final model comparisons.

---

# LASSO performance

| Model | Features | Median AUC | Median accuracy | Non-convergent runs |
|---|---:|---:|---:|---:|
| **Gene** | 3 | **0.990** | **0.941** | 9% |
| **Reaction** | 7 | 0.948 | 0.882 | **0%** |
| **Combined** | 10 | 0.933 | 0.912 | **83%** |

The **gene model achieved the highest predictive performance**.

The reaction model had somewhat lower AUC and accuracy but was the most computationally stable.

The combined model did **not** improve prediction.

Instead, it failed to produce stable coefficients in 83% of the repeated splits.

The third hypothesis was therefore rejected.

More features did not produce a better model.

---

# Important LASSO limitation: information leakage

The reported LASSO performance should be interpreted cautiously.

Candidate features were initially filtered using adjusted paired-Wilcoxon p-values calculated from the complete set of 52 patients **before** the train/test split.

Therefore, the test samples were not completely independent from the feature-selection process.

This represents a potential source of **data leakage** and may lead to optimistic performance estimates.

A stronger implementation would perform all feature filtering independently inside each training split or nested cross-validation workflow.

For this reason, the LASSO models should be treated as:

**exploratory feature-selection and classification models**

rather than externally validated diagnostic models.

---

# Main conclusions

The project produced several main conclusions.

### Transcriptomics and fluxes describe different layers of metabolism

The highest-ranked metabolic subsystems differed substantially between TPM and metabolic-flux analyses.

### Statistical significance does not equal strong classification

Thousands of genes and hundreds of reactions showed statistically significant changes, while substantially fewer achieved high AUC values.

### Integrating analyses reduced the candidate space

Combining DESeq2, TPM and metabolic-flux results reduced thousands of initial features to:

**5 genes associated with 4 reactions.**

### PDAC showed strong patient-to-patient heterogeneity

Both transcriptomic and flux-level analyses showed substantial variation between individual patients.

### More complex models were not necessarily better

The combined LASSO model performed worse and was considerably less stable than the gene-only model.

---

# Limitations

Several limitations should be considered when interpreting the results.

## Limited number of paired patients

Only **52 patients** contained both tumour and normal samples suitable for paired analyses.

This limits statistical power and model validation.

---

## No independent validation dataset

The LASSO models were evaluated using repeated splits of the same underlying cohort.

They were not validated on a completely independent cohort.

---

## Feature-selection leakage

Initial statistical filtering was performed before splitting the paired dataset into training and test sets.

This may inflate estimates of model performance.

---

## Bulk RNA-seq

Bulk RNA-seq represents average expression across all cell populations in a sample.

PDAC tumours contain heterogeneous mixtures of:

- tumour cells;
- stromal cells;
- immune cells;
- other surrounding tissue.

Without cell-composition information, some observed expression differences cannot be assigned specifically to tumour cells.

---

## Limited clinical metadata

The analysis did not include detailed information such as:

- tumour stage;
- patient age;
- sex;
- clinical history;
- tumour morphology;
- cell composition.

These variables could help explain some of the observed heterogeneity.

---

## Model-estimated fluxes

Metabolic fluxes were estimated computationally.

They are not direct measurements of enzyme activity or metabolite concentrations.

---

## Human1 version

The analysis relied on Human1 v1.19.0.

Changes to reaction and gene annotations in later versions can alter biological interpretation.

---

## Choice of distance metric

Metabolic subsystem change was evaluated using normalized Euclidean distance.

Other distance metrics may rank subsystem changes differently.

---

# Future work

Potential extensions include:

- validation using an independent dataset;
- performing all feature selection within training data;
- external validation of LASSO models;
- experimental validation of candidate enzymes and reactions;
- comparison with single-cell RNA-seq;
- integration of clinical and patient metadata;
- integration of cell-composition information;
- repeating the analysis using the newest Human1 version;
- comparing Euclidean distance with alternative metrics such as Manhattan distance;
- experimentally measuring enzyme activity or metabolic flux;
- investigating the identified candidate features in larger PDAC cohorts.

---

# Tools

The analyses were primarily performed in **R**.

### Environment

- R 4.5.2
- RStudio 2026.1.1.403
- Human1 v1.19.0

### Main R packages

- `DESeq2`
- `apeglm`
- `dplyr`
- `tidyr`
- `stringr`
- `readr`
- `readxl`
- `biomaRt`
- `pROC`
- `glmnet`

These were used for:

- data cleaning;
- gene annotation;
- differential-expression analysis;
- statistical testing;
- ROC analysis;
- metabolic-data integration;
- feature selection;
- logistic regression;
- visualization.

---

# Repository structure

```text
PDAC_thesis/
│
├── analyses_pipelines/
│   └── Analysis scripts and intermediate research code
│
├── progress_presentation/
│   └── Presentations documenting the development of the project
│
├── Jakob_Urh_Thesis.pdf
│   └── Final bachelor's thesis
│
├── Thesis_presentation.pdf
│   └── Final thesis defense presentation
│
└── README.md
```

---

# Why the intermediate material is included

This repository intentionally contains more than the final polished thesis.

The folders:

```text
analyses_pipelines/
progress_presentation/
```

preserve parts of the actual research process.

They may therefore contain:

- exploratory code;
- intermediate results;
- preliminary figures;
- approaches that were later changed;
- analyses that were eventually abandoned;
- figures or interpretations that were superseded by the final thesis.

I chose to preserve these files because research is rarely a straight path from a question to a finished result.

The intermediate material documents:

**what I tried, what changed, what failed and how the final analysis developed.**

The authoritative final results are those reported in:

- `Jakob_Urh_Thesis.pdf`
- `Thesis_presentation.pdf`

Intermediate presentations should therefore be interpreted as snapshots of the project at different stages rather than final conclusions.

---

# Research code

The code in `analyses_pipelines/` is intentionally preserved relatively close to the versions used during the thesis.

It is not presented as a polished software package.

Some scripts may therefore contain:

- exploratory sections;
- repeated operations;
- intermediate visualizations;
- research-specific assumptions;
- code written while the analytical strategy was still evolving.

I decided not to rewrite the entire project after completion simply to make the research process appear cleaner than it actually was.

The repository is intended to show both:

**the final analysis**

and

**how I learned to perform it.**

---

# Reproducibility

This repository should not currently be interpreted as a one-command reproducible pipeline.

The original transcriptomic datasets and metabolic-flux input files are not redistributed here.

Reproducing the complete analysis requires:

- the corresponding CPTAC-3 data;
- TCGA-PAAD data;
- metabolic-flux estimates used in the thesis;
- Human1 v1.19.0;
- the required R packages;
- adaptation of local file paths where necessary.

Exact methodological details and data sources are documented in the full thesis.

---

# What I learned

This project was considerably different from my earlier independent data-analysis projects.

## Working with biological data requires context

A statistically strong feature is not automatically biologically meaningful.

Genes, proteins, reactions and metabolic pathways represent different biological layers.

Connecting them requires understanding both the data and the biological model.

---

## Statistical significance and usefulness are different

Many genes and reactions were statistically different between tumour and normal samples.

Far fewer showed strong discriminatory performance.

This made the difference between:

**p-value**

and

**effect / predictive usefulness**

much more concrete for me.

---

## Integration creates new problems

Combining transcriptomic and metabolic results was not simply a matter of joining two tables.

It required dealing with:

- different identifier systems;
- many-to-many gene-reaction relationships;
- missing mappings;
- multiple metabolic subsystems;
- different scales;
- different biological interpretations.

---

## Models are assumptions, not reality

Human1 and estimated metabolic fluxes provide useful abstractions of cellular metabolism.

They are still models.

The SLC16A3–MAR06022 mapping changing between Human1 versions was a particularly useful reminder that computational conclusions can depend on the assumptions encoded in the model.

---

## Patient heterogeneity matters

The paired design made it possible to see how differently individual patients changed between normal and tumour tissue.

Average results can hide substantial variation.

---

## Validation matters

The LASSO analysis produced very high AUC values.

At the same time, examining the methodology revealed limitations including:

- a small dataset;
- feature-selection leakage;
- lack of external validation;
- instability of the combined model.

This changed how I interpret apparently impressive machine-learning metrics.

A high AUC is not enough by itself.

---

## Rejecting a hypothesis is still a result

Two of the three original hypotheses were not supported.

That did not make the project unsuccessful.

It showed that:

- transcript abundance and metabolic flux do not necessarily agree;
- combining more biological layers does not automatically improve prediction;
- biological systems are often more complicated than the initial hypothesis.

That is one of the most important lessons I took from the project.

---

# Research context

This project was completed as my bachelor's thesis in Biochemistry at the **University of Ljubljana, Faculty of Chemistry and Chemical Technology**.

The work was performed under the mentorship of:

**Prof. Dr. Miha Moškon**

with working mentorship from:

**Assistant Ana Halužan Vasle**

The metabolic-flux input data originated from previous work by **Filip Petrovič** and were subsequently analysed and integrated with transcriptomic data in this thesis.

For the full scientific background, references, methodology and discussion, see:

**[Jakob_Urh_Thesis.pdf](Jakob_Urh_Thesis.pdf)**

---

# Final note

This repository represents both the final result of my bachelor's research and the process that produced it.

The polished thesis shows where the project ended.

The code and progress presentations show how it got there.
