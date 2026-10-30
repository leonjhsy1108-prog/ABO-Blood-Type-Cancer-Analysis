# Association Between ABO Blood Type and Cancer: A Statistical Analysis

## Abstract

This study examines the association between ABO blood type and cancer using Pearson chi-square tests of independence on published phenotype counts for 2,640 cancer patients and 2,640 controls. Two levels are analyzed: the observed **blood phenotype** distribution (A, B, AB, O) and the **allele frequencies** (A, B, O) inferred from those phenotypes under Hardy–Weinberg equilibrium. The phenotype distribution of the pooled cancer group differs strongly from controls (χ² = 233.7, df = 3, p < 10⁻¹⁶), driven mainly by an excess of type A and a deficit of type O among cancer patients. Nine of twelve individual cancer types show significant phenotype and allele-level differences from controls; gastrointestinal tract, gall bladder, and (at the allele level) lymphoma do not. The analysis was originally conducted in R and has been fully reimplemented in Python; the original R-based report is included as `ABO Blood Type Analysis.pdf`.

## 1. Introduction

The ABO blood group system, determined by antigens on the surface of red blood cells, is best known for its role in transfusion compatibility. Epidemiological studies have also reported associations between specific blood groups and malignancies such as gastric, pancreatic, and ovarian cancers, raising questions about the biological mechanisms linking blood group antigens and tumor development.

This study formally assesses whether:

1. Blood phenotype distribution differs between cancer patients and controls.
2. Inferred ABO allele frequencies differ between cancer patients and controls.
3. Individual cancer types show distinct blood type associations.

## 2. Data

The counts come from a published case–control study of ABO blood groups in cancer patients, available from the *Internet Journal of Pathology* (Internet Scientific Publications): <https://ispub.com/IJPA/13/1/5982>. The data are aggregate counts, not individual records, and are embedded directly in the analysis script.

| Group | N | A | B | AB | O |
| --- | ---: | ---: | ---: | ---: | ---: |
| Controls | 2,640 | 490 (18.6%) | 1,070 (40.5%) | 300 (11.4%) | 780 (29.5%) |
| All cancers | 2,640 | 903 (34.2%) | 1,068 (40.5%) | 246 (9.3%) | 423 (16.0%) |

Twelve cancer types are reported separately: oral cavity (n = 560), breast (462), cervix (406), ovary (294), gastrointestinal tract (210), bone (189), lungs (168), central nervous system (112), lymphoma (84), gall bladder (70), urinary bladder (50), and kidney (35).

## 3. Methods

### 3.1 Hypotheses

For every comparison:

- **H₀:** blood type and cancer status are independent.
- **H₁:** blood type and cancer status are associated.

All tests use α = 0.05.

### 3.2 Phenotype analysis

A Pearson chi-square test of independence (without continuity correction) is applied to a 2 × 4 contingency table of controls versus the pooled cancer group, and then to a 2 × 4 table of controls versus each individual cancer type.

### 3.3 Allele-frequency analysis

Phenotype counts do not reveal genotypes directly (type A may be AA or AO, type B may be BB or BO). Allele frequencies (p<sub>A</sub>, p<sub>B</sub>, p<sub>O</sub>) are therefore estimated by **maximum likelihood under Hardy–Weinberg equilibrium**, using a grid search with step 0.01 over the multinomial likelihood of the four observed phenotype counts. Expected allele counts (2N × frequency) are then compared between controls and cancer patients in a 2 × 3 chi-square test, overall and for each cancer type.

Estimated allele frequencies:

| Group | p<sub>A</sub> | p<sub>B</sub> | p<sub>O</sub> |
| --- | ---: | ---: | ---: |
| Controls | 0.16 | 0.31 | 0.53 |
| All cancers | 0.26 | 0.30 | 0.44 |

## 4. Results

### 4.1 Phenotype: pooled cancer group versus controls

χ² = 233.73, df = 3, p < 10⁻¹⁶. The phenotype distribution differs significantly. The direction is visible in the data table above: type A is nearly twice as frequent among cancer patients as among controls (34.2% vs 18.6%), type O is about half as frequent (16.0% vs 29.5%), and types B and AB are similar.

### 4.2 Phenotype: individual cancer types versus controls

| Cancer type | n | χ² (df = 3) | p-value | Significant at 0.05 |
| --- | ---: | ---: | ---: | :---: |
| Breast | 462 | 131.90 | < 10⁻¹⁶ | Yes |
| Lungs | 168 | 114.89 | < 10⁻¹⁶ | Yes |
| Ovary | 294 | 90.65 | < 10⁻¹⁶ | Yes |
| Oral cavity | 560 | 80.86 | < 10⁻¹⁶ | Yes |
| Cervix | 406 | 71.19 | 2.4 × 10⁻¹⁵ | Yes |
| Urinary bladder | 50 | 52.42 | 2.4 × 10⁻¹¹ | Yes |
| Kidney | 35 | 39.02 | 1.7 × 10⁻⁸ | Yes |
| Central nervous system | 112 | 38.42 | 2.3 × 10⁻⁸ | Yes |
| Bone | 189 | 26.73 | 6.7 × 10⁻⁶ | Yes |
| Lymphoma | 84 | 23.58 | 3.1 × 10⁻⁵ | Yes |
| Gall bladder | 70 | 4.16 | 0.245 | No |
| Gastrointestinal tract | 210 | 3.23 | 0.358 | No |

### 4.3 Allele frequencies: pooled cancer group versus controls

χ² = 170.62, df = 2, p < 10⁻¹⁶. The inferred allele distribution differs significantly, with a higher A-allele frequency and lower O-allele frequency in the cancer group.

### 4.4 Allele frequencies: individual cancer types versus controls

| Cancer type | n | χ² (df = 2) | p-value | Significant at 0.05 |
| --- | ---: | ---: | ---: | :---: |
| Lungs | 168 | 144.04 | < 10⁻¹⁶ | Yes |
| Breast | 462 | 102.11 | < 10⁻¹⁶ | Yes |
| Cervix | 406 | 96.37 | < 10⁻¹⁶ | Yes |
| Ovary | 294 | 56.61 | 5.1 × 10⁻¹³ | Yes |
| Oral cavity | 560 | 53.07 | 3.0 × 10⁻¹² | Yes |
| Kidney | 35 | 36.39 | 1.3 × 10⁻⁸ | Yes |
| Urinary bladder | 50 | 31.72 | 1.3 × 10⁻⁷ | Yes |
| Central nervous system | 112 | 24.41 | 5.0 × 10⁻⁶ | Yes |
| Bone | 189 | 15.68 | 3.9 × 10⁻⁴ | Yes |
| Lymphoma | 84 | 4.11 | 0.128 | No |
| Gastrointestinal tract | 210 | 1.27 | 0.531 | No |
| Gall bladder | 70 | 0.73 | 0.694 | No |

The two levels of analysis agree closely. Gastrointestinal tract and gall bladder cancers are non-significant under both; lymphoma is significant at the phenotype level but not at the allele level, which is consistent with its small sample (n = 84) and with the allele test having one fewer degree of freedom.

## 5. Interpretation and limitations

The results indicate robust statistical associations between ABO blood type and cancer presence in this dataset, concentrated in an over-representation of type A (and the A allele) and an under-representation of type O among cancer patients. However:

- **Association, not causation.** The chi-square test does not account for confounders such as age, sex, ethnicity, or geography, and the control group's comparability to the cancer group is not known from aggregate counts.
- **No effect sizes.** Odds ratios and confidence intervals were not estimated; large χ² values partly reflect the large sample.
- **Multiple testing.** 26 tests were run without correction. Under a Bonferroni-corrected threshold of 0.05 / 26 ≈ 0.002, every result significant at 0.05 remains significant (the largest significant p-value is 3.9 × 10⁻⁴), so the conclusions are unchanged, but the correction should be applied formally in any extension.
- **Inferred, not observed, alleles.** The allele-level test treats 2N maximum-likelihood allele counts as if they were observed, which assumes Hardy–Weinberg equilibrium and ignores estimation uncertainty. This inflates the allele-level χ² statistics, so they should be read as a consistency check on the phenotype analysis rather than as independent evidence.
- **Small subgroups.** Kidney (n = 35) has an expected AB count below 5 and urinary bladder (n = 50) has no observed AB cases, conditions under which the chi-square approximation is weak.

## 6. Possible extensions

- Logistic regression or odds ratios with confidence intervals for each blood type.
- Post hoc comparisons to identify which phenotype drives each significant result.
- Formal multiple-testing correction (Bonferroni or Benjamini–Hochberg).
- Fisher's exact test for the small subgroups.

## 7. Reproducing the analysis

### Requirements

Python 3.8+ with NumPy, SciPy, and pandas (pandas is optional; it only formats the output tables).

```bash
pip install -r requirements.txt
```

### Run

```bash
python abo_bloodtype_cancer_analysis.py
```

The script prints the phenotype and allele-level tests for the pooled group and for each cancer type, with χ², degrees of freedom, p-value, and the decision at α = 0.05. The data are hard-coded in the script, so no download is needed.

### Repository layout

```text
ABO-Blood-Type-Cancer-Analysis/
├── README.md
├── requirements.txt
├── abo_bloodtype_cancer_analysis.py   # full Python reimplementation
└── ABO Blood Type Analysis.pdf        # original R-based report
```

Note that the allele-level χ² values here differ slightly from those in the original R report because the Python implementation estimates allele frequencies with a 0.01 grid search rather than a continuous optimizer. The conclusions are the same.

## Disclaimer

This project is a statistical and methodological exercise on published aggregate data. Its findings are not medical advice and should not be used for clinical decisions.
