"""Fill A13 placeholders in the working manuscript. Does not edit the read-only draft."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "manuscript/bmc/A_main_v3_work.md"

P12 = (
    "Top-1 accuracy and native-truth top-1 accuracy, in that order, were "
    "0.689 and 0.803 for the baseline, 0.680 and 0.803 for the scaled pure model with C = 0.03, "
    "0.666 and 0.803 for the scaled pure model with C = 0.15, 0.728 and 0.789 for the unscaled mixture model, "
    "and 0.728 and 0.789 for HostMix-TOO in MET500 (evaluation n = 437; native-truth n = 76); "
    "0.682 and 0.868, 0.705 and 0.857, 0.670 and 0.868, 0.777 and 0.846, and 0.756 and 0.802 in POG570 "
    "(n = 512; native-truth n = 91); and 0.527 and 0.473, 0.494 and 0.429, 0.491 and 0.454, "
    "0.706 and 0.727, and 0.698 and 0.730 in the auxiliary RNA-seq set (n = 729; native-truth n = 315)."
)

ACK = (
    "This work would not be possible without the participation of our patients and families, "
    "the POG team, and the generous support of the BC Cancer Foundation and Genome British Columbia "
    "(project B20POG). We also acknowledge contributions towards equipment and infrastructure from "
    "Genome Canada and Genome BC (projects 202SEQ, 212SEQ, 12002), Canada Foundation for Innovation "
    "(projects 20070, 30981, 30198, 33408) and the BC Knowledge Development Fund."
)


def main() -> None:
    text = WORK.read_text()
    replacements = [
        (
            "{{P1: Toil value scale and back-transformation, from the matrix-building code}}",
            "by max(2^x − 0.001, 0). The downloaded Toil matrix is on the log2(TPM + 0.001) scale (negative values mark that encoding)",
        ),
        (
            "The gene universe G consisted of the {{P3: 17552}} protein-coding genes of the HGNC complete set [@Seal2023] {{P3: exact definition, e.g. \"that were present in the TCGA/GTEx, MET500 and POG570 matrices\"}} after Ensembl identifiers were mapped to symbols with GENCODE v23 annotation.",
            "The gene universe G consisted of the 17552 HGNC protein-coding symbols [@Seal2023] present in the TCGA/GTEx Toil matrix, the MET500 matrix and the POG570 matrix after Ensembl identifiers were mapped to symbols with GENCODE v23 annotation.",
        ),
        (
            "the 5000 genes with the highest variance {{P2: of z or of log2 TPM}} across TCGA-train",
            "the 5000 genes with the highest variance of these rank-normal scores across TCGA-train",
        ),
        ("{{P4: 26}}", "26"),
        ("{{P6a: in-pool mean, baseline, ρ = 0.2}}", "0.431"),
        ("{{P6b: in-pool mean, HostMix-TOO, ρ = 0.2}}", "0.007"),
        (
            "{{P22: 0.933 and 0.949, n = 373}}",
            "0.933 for the baseline and 0.949 for HostMix-TOO (n = 373)",
        ),
        ("{{P8a: gated model native-truth top-1, expected 0.473}}", "0.473"),
        ("{{P8b: gated model evaluation top-1, expected 0.550}}", "0.550"),
        ("{{P8c: 22-tissue model pool-out host-attraction rate, expected 0.000}}", "0.000"),
        ("{{P8a}}", "0.473"),
        ("{{P8b}}", "0.550"),
        ("{{P8c}}", "0.000"),
        (
            "{{P8d: verify from the stored gate decisions; delete this clause if it cannot be shown}}",
            "(post hoc count of the recorded gate decisions: the baseline was chosen for 552 of 729 evaluation samples and 303 of 315 native-truth samples)",
        ),
        (
            '{{P9: "most from the liver" only if site counts confirm}}',
            "most from the liver (50 of 56)",
        ),
        ("{{P10: 76}}", "76"),
        ("{{P11a: baseline, n = 536}}", "0.724"),
        ("{{P11b: HostMix-TOO, n = 536}}", "0.455"),
        ("{{P14: expected 0.000}}", "0.000"),
        (
            "{{P12: one sentence on top-1 and native-truth accuracy of the five ablation models, with values, or delete this placeholder if those metrics were not computed}}",
            P12,
        ),
        ("{{P15: 145 of 147}}", "145 of 147"),
        ("{{P15: cohort}}", "GSE41258"),
        ("{{P7: samples in the 22 tissues}}", "1874"),
        (
            "{{P16: totals for the one-per-patient sets, e.g. 536 / 207 / 341}}",
            "536 / 207 / 341",
        ),
        (
            "{{P18: SCOPE commit or \"not recorded\"}}",
            "fd9db4cd20122c8c7b675400c0df2a4c5bebb257",
        ),
        ("{{ACK-POG: POG570 acknowledgement as required by the POG570 README; see A 지시서 13 §6.2}}", ACK),
        ("[{{P17: data paper of GSE14018}}", "[@Zhang2009"),
        ("[{{P17: data paper of GSE71729}}", "[@Moffitt2015"),
        ("[{{P17: data paper of GSE74685}}", "[@Haider2015"),
        (
            "| GSE41258 | Descriptive, microarray | GPL96 | Colorectal cancer | {{P16}} | [@Sheffer2009] |",
            "| GSE41258 | Descriptive, microarray | GPL96 | Colorectal cancer | 236 / 58 / 183 | [@Sheffer2009] |",
        ),
        (
            "| GSE14018 | Descriptive, microarray | GPL96 | Breast cancer | {{P16}} | [{{P17}}] |",
            "| GSE14018 | Descriptive, microarray | GPL96 | Breast cancer | 36 / 36 / 0 | [@Zhang2009] |",
        ),
        (
            "| GSE71729 | Descriptive, microarray | GPL20769 | Pancreatic cancer | {{P16}} | [{{P17}}] |",
            "| GSE71729 | Descriptive, microarray | GPL20769 | Pancreatic cancer | 206 / 56 / 145 | [@Moffitt2015] |",
        ),
        (
            "| GSE74685 | Descriptive, microarray | GPL15659 | Prostate cancer | {{P16}} | [{{P17}}] |",
            "| GSE74685 | Descriptive, microarray | GPL15659 | Prostate cancer | 38 / 49 / 0 | [@Haider2015] |",
        ),
        (
            "| FHCRC | Descriptive, microarray | Agilent | Prostate cancer | {{P16}} | [@Kumar2016] |",
            "| FHCRC | Descriptive, microarray | Agilent | Prostate cancer | 20 / 8 / 13 | [@Kumar2016] |",
        ),
        (
            "For the two-channel array with a common reference, the sample channel was used.",
            "For the two-channel arrays GSE71729 and GSE74685, the sample channel was used and the reference channel was not.",
        ),
        (
            "The expression file was hashed at download and could not be read by the analysis code until the POG570 analysis plan had been committed.",
            "The expression file was hashed at download. The analysis code returns that matrix only when called with an explicit unlock, and that call was made only after the POG570 analysis plan had been committed.",
        ),
        (
            "we regenerated the HostMix-TOO training mixtures from the stored records and refitted HostMix-TOO",
            "we regenerated the HostMix-TOO training mixtures from the mixture records and refitted HostMix-TOO",
        ),
        (
            "Using stored predictions and frozen models only",
            "Using the confirmatory predictions and the locked models only",
        ),
    ]
    for old, new in replacements:
        count = text.count(old)
        if count != 1:
            raise SystemExit(f"expected 1 occurrence, found {count}: {old[:80]}")
        text = text.replace(old, new)
    text = text.split("-->", 1)[-1].lstrip("\n") if text.startswith("<!--") else text
    leftover = [line for line in text.splitlines() if "{{" in line]
    if leftover:
        raise SystemExit("placeholders remain:\n" + "\n".join(leftover[:12]))
    WORK.write_text(text)
    print("filled", WORK)


if __name__ == "__main__":
    main()
