
# Association Between ABO Blood Type and Cancer: A Statistical Analysis 
## Abstract

This study examines the association between ABO blood type and various cancer types using chi-square tests of independence. Both blood phenotype (A, B, AB, O) and genotype distributions were analyzed in relation to cancer presence. Statistical evidence indicates significant associations between blood type and most cancer types examined. The analysis was originally conducted in R and has been fully reimplemented in Python to ensure reproducibility and transparency. The original R-based report is available in the file ABO Blood Type Analysis.pdf included in this repository.

---

## 1. Introduction

The ABO blood group system, determined by antigens present on the surface of red blood cells, has long been recognized for its importance in transfusion medicine. Beyond its clinical role in compatibility, epidemiological studies have suggested that blood type may influence susceptibility to certain diseases, including cancer.

Prior research has reported associations between specific blood groups and malignancies such as gastric, pancreatic, and ovarian cancers. These findings raise important questions regarding the biological mechanisms linking blood group antigens and tumor development.

The purpose of this study is to formally assess whether:

1. Blood phenotype distribution differs between cancer patients and controls.
2. Blood genotype distribution differs between cancer patients and controls.
3. Individual cancer types demonstrate distinct blood type associations.

Data source:
https://ispub.com/IJPA/13/1/5982

---

## 2. Statistical Methodology

### 2.1 Hypothesis Testing Framework

For both phenotype and genotype analyses, the following hypotheses were tested:

- **Null Hypothesis (H₀):**  
  Blood type and cancer presence are independent.

- **Alternative Hypothesis (H₁):**  
  Blood type and cancer presence are associated.

All tests were conducted at a significance level of α = 0.05.

### 2.2 Statistical Procedure

The primary method used was the Pearson Chi-Square Test of Independence, applied to contingency tables constructed from:

- Control group blood type counts  
- Cancer group blood type counts  

For genotype analysis, genotype frequencies were approximated via maximum likelihood estimation derived from phenotype counts.

---

## 3. Phenotype Analysis

### 3.1 Overall Cancer vs. Control Comparison

A 2×4 contingency table was constructed comparing phenotype distributions between controls and the aggregated cancer group.

Results:

- χ² = 233.73  
- df = 3  
- p-value < 2.2 × 10⁻¹⁶  

Conclusion: Blood phenotype distribution differs significantly between cancer patients and controls.

---

### 3.2 Individual Cancer Type Comparisons (Phenotype)

Statistically significant associations were observed for:

- Oral cavity
- Lungs
- Breast
- Cervix
- Ovary
- Lymphoma
- Bone
- Kidney
- Urinary bladder
- Central nervous system

Non-significant associations:

- Gastrointestinal tract
- Gall bladder

---

## 4. Genotype Analysis

Genotype frequencies were approximated using maximum likelihood estimation based on observed phenotype counts.

### 4.1 Overall Genotype Comparison

- χ² = 165.88  
- df = 2  
- p-value < 2.2 × 10⁻¹⁶  

Conclusion: Blood genotype distribution differs significantly between cancer patients and controls.

---

### 4.2 Individual Cancer Type Comparisons (Genotype)

All cancer types exhibited statistically significant differences in genotype distribution relative to controls.

---

## 5. Interpretation and Limitations

The results indicate robust statistical associations between ABO blood type and cancer presence.

However:

- The chi-square test identifies association, not causation.
- Effect size and relative risk were not estimated.
- Specific blood types contributing to increased susceptibility were not isolated.
- Multiple testing adjustments were not applied.

Further modeling would be necessary to quantify risk and determine clinical relevance.

---

## 6. Potential Extensions

Future analyses could incorporate:

- Logistic regression modeling  
- Odds ratio estimation  
- Post hoc pairwise comparisons  
- Multiple testing correction  
- Effect size estimation  

---

## 7. Implementation Notes

This repository contains a complete Python reimplementation of the original R-based analysis. The code reproduces:

- Contingency table construction  
- Chi-square tests of independence  
- Maximum likelihood estimation for genotype frequencies  
- Individual cancer-type comparisons  

---

## 8. Reproducibility

Install dependencies:

    pip install -r requirements.txt

Run the analysis:

    python abo_bloodtype_cancer_analysis.py

---

## 9. Summary of Findings

- Strong statistical evidence of association between blood phenotype and cancer.
- Strong statistical evidence of association between blood genotype and cancer.
- Nearly all cancer types show significant associations.
- Gastrointestinal tract and gall bladder cancers were not significant in phenotype analysis.
- Genotype analysis showed consistent significance across all cancer types.

---

## Disclaimer

This project is intended for statistical and methodological exploration.  
Findings should not be interpreted as medical advice or clinical guidance.
