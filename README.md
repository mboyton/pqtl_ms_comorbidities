# Plasma proteome Mendelian randomisation analysis in multiple sclerosis

This repository contains analysis code used to perform **proteome-wide Mendelian randomisation (MR)** and **colocalisation analyses** to identify plasma proteins with potential causal effects on **multiple sclerosis (MS)** risk.

The analysis integrates **plasma protein quantitative trait loci (pQTL)** data with **genome-wide association study (GWAS)** summary statistics to prioritise candidate therapeutic targets and evaluate potential cross-disease effects across common MS comorbidities.

---

# Repository structure

```
.
├── pqtl_ms_and_comorbidities_mendelian_randomisation.R
├── colocalisation_decode_pqtl_ms.R
└── README.md
```

---

# Analysis overview

The analysis consists of three main stages.

## 1. Proteome-wide Mendelian randomisation

Proteome-wide two-sample Mendelian randomisation is performed using **cis-pQTLs** as instruments to estimate the causal effect of circulating plasma proteins on MS risk.

Proteins with **FDR-adjusted p < 0.05** are prioritised for downstream analyses.

Script:

```
pqtl_ms_and_comorbidities_mendelian_randomisation.R
```

---

## 2. Cross-phenotype Mendelian randomisation

Proteins identified in the primary MS MR screen are tested against six clinically relevant comorbid conditions.

Comorbid outcomes analysed:

- Major depression
- Generalised anxiety disorder
- Hypertension
- Hypercholesterolaemia
- Asthma
- Chronic obstructive pulmonary disease (COPD)

Outcome summary statistics are retrieved from **OpenGWAS** using the `TwoSampleMR` package.

---

## 3. Colocalisation analysis

To distinguish true causal relationships from confounding due to linkage disequilibrium, **colocalisation analysis** is performed using the `coloc` R package.

Approximate Bayes Factor (ABF) colocalisation estimates the posterior probability that the **same causal variant** drives both the pQTL and MS GWAS signal.

Posterior probabilities are estimated for five hypotheses:

| Hypothesis | Interpretation |
|------------|---------------|
| H0 | No association with either trait |
| H1 | Association with protein only |
| H2 | Association with MS only |
| H3 | Both traits associated but with distinct causal variants |
| H4 | Shared causal variant |

Strong evidence of colocalisation in our analyses is defined as:

```
PP.H4 > 0.8
```

Script:

```
colocalisation_decode_pqtl_ms.R
```

---

# Data sources

## Plasma proteome

**deCODE genetics plasma pQTL dataset**

Ferkingstad E, Sulem P, Atlason BA, et al. Large-scale integration of the plasma proteome with genetics and disease. Nat Genet. Dec 2021;53(12):1712-1721. doi:10.1038/s41588-021-00978-w

---

## Multiple sclerosis GWAS

International Multiple Sclerosis Genetics Consortium (IMSGC)

Patsopoulos NA, Baranzini SE, et al. Multiple sclerosis genomic map implicates peripheral immune cells and microglia in susceptibility. Science. 2019;365(6460):eaav7188. doi:10.1126/science.aav7188

OpenGWAS dataset:

```
ieu-b-18
```

---

## Comorbidity GWAS

OpenGWAS datasets used for cross-trait MR analyses:

| Trait | OpenGWAS ID |
|------|-------------|
| Major depression | ieu-b-102 |
| Generalised anxiety disorder | finn-b-F5_GAD |
| Hypertension | ukb-b-12493 |
| Hypercholesterolaemia | finn-b-E4_HYPERCHOL |
| Asthma | ebi-a-GCST90014325 |
| COPD | finn-b-J10_COPD |

---

# Software

All analyses for the accompanying manuscript were conducted using **R (≥ 4.3)** on the University of Bristol's High Performance Computing cluster.

Required R packages:

```
TwoSampleMR
coloc
dplyr
ggplot2
patchwork
stringr
readxl
SciViews
```

---

# Running the analysis

File paths in the scripts are provided as **placeholders** and should be updated to point to local data files.

---
