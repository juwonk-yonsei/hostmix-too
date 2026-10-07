#!/usr/bin/env python3
"""Build the review-pack contact sheet and the release data manifest."""
from __future__ import annotations

import csv
import shutil
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "manuscript" / "bmc" / "review_pack"

PAGES = [
    ("manuscript/bmc/figures/fig1_simulation.png", "Figure 1. Stored simulation of host-pull rate.\n\nEach point is the unweighted mean of the stored host-pull rate across tissue-site rows at that mixing weight. In-pool sites are those whose native organs sit in the training host pool. Pool-out sites are the twelve sites named in the text. The panel is a description of the stored simulation table, not a preregistered patient-level test."),
    ("manuscript/bmc/figures/fig2_prereg_forest.png", "Figure 2. Preregistered SA-Z minus BASE-Z differences.\n\nPoints are the stored differences. Intervals are the stored patient-bootstrap percentile intervals. POG570 contrasts are H1, H2, and H3. Auxiliary contrasts are AH1–AH3 and the secondary contrasts AS1–AS3."),
    ("manuscript/bmc/figures/fig3_cohort_forest.png", "Figure 3. Cohort meta-analysis of the primary contrasts.\n\nThe vertical line is the stored DerSimonian–Laird estimate. Cohort intervals are the stored effect plus or minus 1.96 times the square root of the stored variance. The seven-cohort RNA-seq set includes GSE209998, which was not part of the six-cohort preregistered family."),
    ("manuscript/bmc/figures/fig4_method_heatmap.png", "Figure 4. Locked development comparison on MET500.\n\nThe heatmap is the stored comparison of representations and model families on the development cohort. It is context for why SA-Z was locked. It is not a second confirmation. Blank cells are stored as missing."),
    ("manuscript/bmc/figures/fig5_external.png", "Figure 5. Public classifiers and post-hoc set splits.\n\nThe top row is the pre-specified descriptive comparison on the common-label set. It is not a hypothesis test. Top-1 denominators are 437, 512, and 729. Host-rate denominators are 361, 378, and 427. The bottom row is post hoc."),
    ("manuscript/bmc/figures/fig6_limits.png", "Figure 6. Post-hoc error shift and missing genes.\n\nPanel A counts patients where BASE-Z was wrong and SA-Z was right. Panel B is the fraction of non-esophagus truths called esophagus. Panel C is the GSE41258 colon and rectum slice, selected evaluation n = 183. Panel D is the change in top-1 after restricting genes to a microarray platform."),
    ("manuscript/bmc/additional_figures/figS1.png", "Additional file 1: Figure S1. IMvigor210 esophagus-call contributions.\n\nMean logistic contribution on IMvigor210 bladder native-truth samples predicted as esophagus. BASE-Z n = 145. SA-Z n = 65. The native-truth denominator is 194. The gene list is not a validated marker panel."),
    ("manuscript/bmc/additional_figures/figS2.png", "Additional file 1: Figure S2. Esophagus calls and ESCA histology centers.\n\nTCGA-train esophagus n = 145, split into adenocarcinoma 74, squamous carcinoma 65, and other 6. Each selected-evaluation sample called esophagus is assigned by the larger Spearman correlation."),
]


def contact_sheet() -> None:
    PACK.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "manuscript/bmc/A_main.md", PACK / "A_main.md")
    shutil.copy2(ROOT / "manuscript/A_supplement.md", PACK / "Additional_file_1.md")
    out = PACK / "figures_contact_sheet.pdf"
    with PdfPages(out) as pdf:
        for relative, legend in PAGES:
            image = plt.imread(ROOT / relative)
            figure = plt.figure(figsize=(8.5, 11))
            axis = figure.add_axes((0.08, 0.28, 0.84, 0.66))
            axis.imshow(image)
            axis.axis("off")
            figure.text(0.08, 0.06, legend, fontsize=9, va="bottom", wrap=True)
            pdf.savefig(figure)
            plt.close(figure)
    print("contact", out.stat().st_size)


def terms(source: str, license_text: str) -> str:
    if source == "POG570":
        return "URL and hash only. The file is not in the release package."
    if license_text:
        return "cBioPortal file. ODbL. See the study LICENSE."
    if source in {"toil", "met500", "xena"}:
        return "Public hub download. Expression matrices are not redistributed in the release package."
    if source == "hgnc":
        return "HGNC public download."
    if source == "geo":
        return "NCBI GEO public download."
    if source == "hpa":
        return "Human Protein Atlas public download."
    if source == "oncfind":
        return "Git clone. Commit recorded in the hash column."
    return "Public download recorded in the project log."


def manifest() -> None:
    rows = []
    log = ROOT / "logs" / "download_log.tsv"
    with log.open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            rows.append({
                "dataset": rec["source"],
                "file": Path(rec["file"]).name,
                "url": rec["url"],
                "bytes": rec["size_bytes"],
                "sha256": rec["sha256"],
                "date": rec["datetime"][:10],
                "terms": terms(rec["source"], ""),
            })
    aux = ROOT / "config" / "aux_downloads_A5.tsv"
    with aux.open() as handle:
        for rec in csv.DictReader(handle, delimiter="\t"):
            local = ROOT / "data/processed/aux/downloads" / f"{rec['cohort']}__{rec['file']}"
            received = ""
            if local.exists():
                received = datetime.fromtimestamp(local.stat().st_mtime).date().isoformat()
            rows.append({
                "dataset": rec["cohort"],
                "file": rec["file"],
                "url": rec["url"],
                "bytes": rec["bytes"],
                "sha256": rec["sha256"],
                "date": received,
                "terms": terms(rec["cohort"], rec["license"]),
            })
    dest = ROOT / "release" / "data_manifest.tsv"
    with dest.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "file", "url", "bytes", "sha256", "date", "terms"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("manifest rows", len(rows))


if __name__ == "__main__":
    contact_sheet()
    manifest()
