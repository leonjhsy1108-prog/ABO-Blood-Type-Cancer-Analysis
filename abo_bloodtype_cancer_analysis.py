from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import chi2_contingency

try:
    import pandas as pd
except Exception:
    pd = None  




PHENOTYPES = ["A", "B", "AB", "O"]
ALLELES = ["A", "B", "O"]

CONTROL_PHENO = {"A": 490, "B": 1070, "AB": 300, "O": 780}
TOTAL_CANCER_PHENO = {"A": 903, "B": 1068, "AB": 246, "O": 423}

CANCER_TYPES = [
    "Oral cavity",
    "Gastrointestinal tract",
    "Gall bladder",
    "Lungs",
    "Breast",
    "Cervix",
    "Ovary",
    "Lymphoma",
    "Bone",
    "Kidney",
    "Urinary bladder",
    "Central nervous system",
]

CANCER_PHENO_MATRIX = np.array(
    [
        [196, 210, 42, 112],
        [49, 84, 21, 56],
        [14, 35, 4, 17],
        [84, 56, 21, 7],
        [196, 140, 28, 98],
        [112, 154, 84, 56],
        [98, 154, 21, 21],
        [28, 42, 2, 12],   
        [42, 105, 14, 28],
        [21, 5, 2, 7],     
        [28, 20, 0, 2],
        [35, 63, 7, 7],
    ],
    dtype=int,
)



def dict_to_row(d: Dict[str, int], labels: List[str]) -> np.ndarray:
    return np.array([int(d[k]) for k in labels], dtype=int)


def chi_square_test(table: np.ndarray) -> Dict[str, float]:

    chi2, p, dof, expected = chi2_contingency(table, correction=False)
    return {"chi2": float(chi2), "dof": int(dof), "p_value": float(p)}


def format_p(p: float) -> str:
    if p < 1e-16:
        return "< 1e-16"
    return f"{p:.4g}"



@dataclass(frozen=True)
class AlleleMLE:
    pA: float
    pB: float
    pO: float
    loglik: float


def allele_mle_from_phenotypes(
    nA: int,
    nB: int,
    nAB: int,
    nO: int,
    tol: float = 0.01,
) -> AlleleMLE:
    best = AlleleMLE(pA=0.0, pB=0.0, pO=1.0, loglik=-np.inf)

    p_grid = np.arange(0.0, 1.0 + 1e-12, tol)

    for pA in p_grid:
        for pB in p_grid:
            if pA + pB >= 1.0:
                continue
            pO = 1.0 - pA - pB

            pA_ph = pA * pA + 2.0 * pA * pO
            pB_ph = pB * pB + 2.0 * pB * pO
            pO_ph = pO * pO
            pAB_ph = 2.0 * pA * pB

            probs = np.array([pA_ph, pB_ph, pAB_ph, pO_ph], dtype=float)
            if np.any(probs <= 0):
                continue

            loglik = (
                nA * np.log(pA_ph)
                + nB * np.log(pB_ph)
                + nAB * np.log(pAB_ph)
                + nO * np.log(pO_ph)
            )

            if loglik > best.loglik:
                best = AlleleMLE(pA=float(pA), pB=float(pB), pO=float(pO), loglik=float(loglik))

    return best


def allele_counts_from_mle(mle: AlleleMLE, N: int) -> np.ndarray:
    counts = 2.0 * N * np.array([mle.pA, mle.pB, mle.pO], dtype=float)
    return np.rint(counts).astype(int)


def run_phenotype_analysis(alpha: float = 0.05) -> None:
    print("=" * 80)
    print("Phenotype analysis (A/B/AB/O) – Pearson chi-square tests")
    print("=" * 80)

    control = dict_to_row(CONTROL_PHENO, PHENOTYPES)
    total_cancer = dict_to_row(TOTAL_CANCER_PHENO, PHENOTYPES)

    table = np.vstack([control, total_cancer])
    res = chi_square_test(table)
    print("\nOverall: Control vs Total Cancer")
    print(f"  X^2 = {res['chi2']:.3f}, df = {res['dof']}, p = {format_p(res['p_value'])}")
    print(f"  Decision at alpha={alpha}: {'REJECT H0' if res['p_value'] < alpha else 'FAIL TO REJECT H0'}")

    print("\nPer-cancer-type: Control vs Cancer Type")
    rows = []
    for i, cancer in enumerate(CANCER_TYPES):
        counts = CANCER_PHENO_MATRIX[i]
        table_i = np.vstack([control, counts])
        r = chi_square_test(table_i)
        rows.append(
            {
                "Cancer type": cancer,
                "chi2": r["chi2"],
                "df": r["dof"],
                "p_value": r["p_value"],
                "p_fmt": format_p(r["p_value"]),
                "sig@0.05": r["p_value"] < alpha,
            }
        )

    if pd is not None:
        df = pd.DataFrame(rows).sort_values("p_value")
        print(df[["Cancer type", "chi2", "df", "p_fmt", "sig@0.05"]].to_string(index=False))
    else:
        for row in rows:
            print(f"  {row['Cancer type']:<26} X^2={row['chi2']:.3f}, df={row['df']}, p={row['p_fmt']}")


def run_genotype_allele_analysis(alpha: float = 0.05, tol: float = 0.01) -> None:
    print("\n" + "=" * 80)
    print("Genotype / allele analysis (A/B/O allele counts inferred from phenotypes)")
    print("=" * 80)
    print(f"Allele MLE grid step (tol) = {tol}")

    control_row = dict_to_row(CONTROL_PHENO, PHENOTYPES)
    total_cancer_row = dict_to_row(TOTAL_CANCER_PHENO, PHENOTYPES)

    def mle_and_counts(row: np.ndarray) -> Tuple[AlleleMLE, np.ndarray, int]:
        nA, nB, nAB, nO = map(int, row.tolist())
        N = nA + nB + nAB + nO
        mle = allele_mle_from_phenotypes(nA, nB, nAB, nO, tol=tol)
        counts = allele_counts_from_mle(mle, N)
        return mle, counts, N

    mle_control, alleles_control, N_control = mle_and_counts(control_row)
    mle_total, alleles_total, N_total = mle_and_counts(total_cancer_row)

    print("\nMLE allele frequencies:")
    print(f"  Control (N={N_control}): pA={mle_control.pA:.3f}, pB={mle_control.pB:.3f}, pO={mle_control.pO:.3f}")
    print(f"  Cancer  (N={N_total}): pA={mle_total.pA:.3f}, pB={mle_total.pB:.3f}, pO={mle_total.pO:.3f}")

    print("\nInferred allele counts (rounded):")
    print(f"  Control: {dict(zip(ALLELES, alleles_control.tolist()))}")
    print(f"  Cancer : {dict(zip(ALLELES, alleles_total.tolist()))}")

    table = np.vstack([alleles_control, alleles_total])
    res = chi_square_test(table)
    print("\nOverall: Control vs Total Cancer (alleles A/B/O)")
    print(f"  X^2 = {res['chi2']:.3f}, df = {res['dof']}, p = {format_p(res['p_value'])}")
    print(f"  Decision at alpha={alpha}: {'REJECT H0' if res['p_value'] < alpha else 'FAIL TO REJECT H0'}")

    print("\nPer-cancer-type: Control vs Cancer Type (alleles A/B/O)")
    rows = []
    for i, cancer in enumerate(CANCER_TYPES):
        row = CANCER_PHENO_MATRIX[i]
        mle_i, alleles_i, N_i = mle_and_counts(row)
        table_i = np.vstack([alleles_control, alleles_i])
        r = chi_square_test(table_i)
        rows.append(
            {
                "Cancer type": cancer,
                "N_type": N_i,
                "pA": mle_i.pA,
                "pB": mle_i.pB,
                "pO": mle_i.pO,
                "chi2": r["chi2"],
                "df": r["dof"],
                "p_value": r["p_value"],
                "p_fmt": format_p(r["p_value"]),
                "sig@0.05": r["p_value"] < alpha,
            }
        )

    if pd is not None:
        df = pd.DataFrame(rows).sort_values("p_value")
        print(df[["Cancer type", "N_type", "chi2", "df", "p_fmt", "sig@0.05"]].to_string(index=False))
    else:
        for row in rows:
            print(f"  {row['Cancer type']:<26} X^2={row['chi2']:.3f}, df={row['df']}, p={row['p_fmt']}")


def main() -> None:
    run_phenotype_analysis(alpha=0.05)
    run_genotype_allele_analysis(alpha=0.05, tol=0.01)


if __name__ == "__main__":
    main()
