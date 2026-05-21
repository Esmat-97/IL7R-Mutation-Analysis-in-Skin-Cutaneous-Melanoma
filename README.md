# IL7R Mutation Analysis in Skin Cutaneous Melanoma (TCGA-SKCM)

## Overview

This project investigates the mutation landscape and potential immunological relevance of the IL7R gene in Skin Cutaneous Melanoma (SKCM) using data from The Cancer Genome Atlas (TCGA) PanCancer Atlas.

The study explores whether IL7R mutations are associated with:

- Tumor Mutation Burden (TMB)
- Immune-related gene expression
- Immune checkpoint expression
- Survival outcomes
- Co-mutation patterns with major melanoma driver genes

---

## Research Question

Does IL7R mutation status define a biologically distinct subset of melanoma characterized by altered immune features, mutation burden, or clinical outcome?

---

## Dataset

Source:

- TCGA Skin Cutaneous Melanoma (SKCM)
- PanCancer Atlas
- cBioPortal
- Cell (2018)

Cohort size:

- 440 mutation-profiled patients
- 444 expression-profiled patients

---

## Analyses Performed

### 1. Mutation Landscape

Identified the most frequently mutated genes in SKCM.

Top mutated genes:

| Gene | Mutation Count |
|--------|--------|
| TTN | 3476 |
| MUC16 | 2274 |
| DNAH5 | 989 |
| PCLO | 771 |
| LRP1B | 578 |

---

### 2. Immune Gene Mutation Analysis

Evaluated mutation frequencies in selected immune-related genes:

- IL7R
- JAK1
- JAK2
- STAT3
- PDCD1
- CD274
- CTLA4
- IFNG
- HLA genes

Key finding:

- IL7R was the most frequently mutated immune-related gene.

| Gene | Mutations |
|--------|--------|
| IL7R | 126 |
| STAT3 | 35 |
| JAK1 | 34 |
| JAK2 | 34 |
| PDCD1 | 30 |

---

### 3. Survival Analysis

Kaplan-Meier and log-rank tests were performed for immune-related genes.

Result:

No significant survival association was detected for IL7R or other analyzed immune genes.

| Gene | Log-rank p-value |
|--------|--------|
| IL7R | 0.856 |
| PDCD1 | 0.464 |
| JAK1 | 0.359 |

---

### 4. IL7R Mutation Spectrum

Variant classifications observed:

| Mutation Type | Count |
|--------|--------|
| Missense Mutation | 56 |
| Silent | 38 |
| 5'UTR | 9 |
| Nonsense Mutation | 7 |
| 3'UTR | 7 |

Missense mutations represented the dominant mutation class.

---

### 5. Co-Mutation Analysis

Investigated overlap between IL7R mutations and major melanoma drivers.

Genes examined:

- BRAF
- NRAS
- TP53

Example:

| Gene | IL7R + Gene |
|--------|--------|
| BRAF | 48 |
| NRAS | 29 |
| TP53 | 26 |

No significant association was detected between IL7R and BRAF mutations.

Fisher's Exact Test:

p = 0.463

---

### 6. Tumor Mutation Burden (TMB)

Compared TMB between IL7R-mutant and wild-type tumors.

Statistical test:

- Mann-Whitney U test

Result:

p = 3.55 × 10⁻¹³

Key finding:

IL7R-mutant tumors exhibit significantly higher mutation burden than wild-type tumors.

---

### 7. Immune Activity Assessment

Evaluated expression of:

- CD8A
- GZMB
- PRF1
- IFNG

Derived scores:

- CD8 score
- Cytolytic score

Results:

No statistically significant increase in immune infiltration signatures was observed in IL7R-mutant tumors.

---

### 8. Chemokine Analysis

Genes:

- CXCL9
- CXCL10

Results:

Borderline trends were observed but did not reach statistical significance.

---

### 9. Immune Checkpoint Analysis

Genes:

- CD274 (PD-L1)
- PDCD1
- CTLA4
- LAG3
- TIGIT
- HAVCR2

Key result:

CD274 (PD-L1) expression showed a modest increase in IL7R-mutant tumors.

| Gene | p-value |
|--------|--------|
| CD274 | 0.047 |
| PDCD1 | 0.165 |
| CTLA4 | 0.390 |

---

### 10. Melanoma Driver Mutation Comparison

Compared mutation counts of:

- BRAF
- NRAS
- NF1

Result:

NF1 displayed a significant difference between IL7R-mutant and wild-type groups.

| Gene | p-value |
|--------|--------|
| NF1 | 0.014 |

---

## Main Findings

- IL7R is the most frequently mutated immune-related gene in TCGA-SKCM.
- IL7R mutations are strongly associated with increased Tumor Mutation Burden.
- IL7R mutation status does not significantly affect overall survival.
- Immune infiltration signatures show limited differences between mutant and wild-type tumors.
- PD-L1 expression is modestly elevated in IL7R-mutant samples.
- NF1 mutation patterns differ between IL7R-mutant and wild-type groups.

---

## Tools

- Python
- Pandas
- NumPy
- Matplotlib
- Seaborn
- SciPy
- Lifelines

---

## Future Directions

Potential extensions include:

- Validation in independent melanoma cohorts
- Analysis of immunotherapy-treated patients
- Functional characterization of recurrent IL7R mutations
- Integration with immune subtype classifications

---

## Author

Mohamed Esmat

Cancer Biology | Immunology | Bioinformatics
