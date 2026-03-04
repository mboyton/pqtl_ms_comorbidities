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

Proteome-wide two-sample Mendelian randomisation is performed using **cis-acting plasma pQTLs** as instruments to estimate the causal effect of circulating proteins on MS risk.

Key steps include:

- formatting pQTL exposure data
- harmonising exposure and outcome alleles
- performing inverse variance weighted (IVW) Mendelian randomisation
- applying false discovery rate (FDR) correction across proteins

Proteins with **FDR-adjusted p < 0.05** are prioritised for downstream analyses.

Script:

```
pqtl_ms_and_comorbidities_mendelian_randomisation.R
```

---

## 2. Cross-phenotype Mendelian randomisation

Proteins identified in the primary MS MR screen are tested against six clinically relevant comorbid conditions to evaluate potential pleiotropic effects.

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

Strong evidence of colocalisation is defined as:

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

Ferkingstad et al. (2021)  
*Large-scale integration of the plasma proteome with genetics and disease*  
Nature Genetics

---

## Multiple sclerosis GWAS

International Multiple Sclerosis Genetics Consortium (IMSGC)

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

Analyses were conducted using **R (≥ 4.3)**.

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

Typical workflow:

```
1. Run pqtl_ms_and_comorbidities_mendelian_randomisation.R
   - performs proteome-wide MR
   - identifies significant proteins
   - performs comorbidity MR analyses

2. Run colocalisation_decode_pqtl_ms.R
   - performs ABF colocalisation for prioritised proteins
```

---

# Reproducibility notes

- Analyses were originally conducted on a **high-performance computing cluster**.
- Scripts are designed to be run sequentially.
- Some datasets (e.g. deCODE pQTL summary statistics) may require controlled access.

---

# Citation

If using this code, please cite the associated manuscript describing this analysis.