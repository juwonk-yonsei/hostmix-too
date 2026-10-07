"""Package release files, supplementary tables, and the history bundle.

Does not rewrite git history and does not include POG570 sample-level files.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from openpyxl import Workbook

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release"
BMC = ROOT / "manuscript" / "bmc"

HISTORY = [
    ("e579bc2d20ce5bae611f18aafe24ebb11f47d742", "A3", "preregistration", "POG570 confirmation plan committed before feature calculation."),
    ("ad5fcb8e017c642904b252deaa5b113dc34b0d63", "A3-A4", "freeze", "POG570 predictions and MET500 conformal freeze."),
    ("6f91d0adb310d3b594d569d14d4b8ba83688ec57", "A4", "confirmation output", "POG570 confirmation output hash. The run was not repeated."),
    ("ba3c074e96e0cfd0124e54a24e08964849e7879e", "A6", "preregistration", "Auxiliary hypothesis preregistration."),
    ("1b5c0f69f117436281ab638e4042340bd9e797e6", "A6", "model freeze", "Stage-6 model, feature, LD, and label freeze."),
    ("84086a1a48f74d0d457ad05244727473c1b8169d", "A7", "code freeze", "Auxiliary confirmation code freeze."),
    ("198007aa", "A7", "output path", "Output-path commit."),
    ("edd32463afe375e87de245beba6969b94ff557fa", "A7", "unlock", "Unlock after the code freeze."),
    ("a7272492", "A7", "deviation fix", "First --real run wrote no hypothesis table because standard_site was empty. This commit is the fix. config/unlock_A7.md was not edited."),
    ("8f788feb74e51961b3e1c60b5dc25165f68cdc6f", "A7", "confirmation output", "Second --real run. This is the preregistered auxiliary result."),
    ("3f2a322e5f14687d423fad40af253bbfb3a9c240", "A7", "addendum", "Chronology of the empty-site failure. The unlock file itself was not edited."),
    ("95dbcf78b8487e66cc3d7a26c2068ee148e4e052", "A8", "preregistration", "Stage-8 comparison plan committed before any cohort was sent to SCOPE or CUP-AI-Dx."),
    ("011b4dbaf33f55053f5b7c92d639255b21a4940b", "A8", "input freeze", "Evaluation-matrix hashes. The 18 GSE50760 risk-only samples are outside this file and were scored later so the host-rate denominator could stay 427."),
    ("fd712e02b17ece1ec7af3c6bfde5501994716ae5", "A7-A8", "label map", "External-tool label map committed before example output was inspected. He et al. 2023 was not run."),
    ("177ca6980a9da42d4fdeac16c807942af76fe5bd", "A8", "tool output hash", "External-classifier output hashes."),
    ("111eee8ff9b29abccf79faadf0a04ea4ae847b5b", "A12", "post hoc plan", "Ablation plan with its reading rule, committed before any ablation model was fitted. Post hoc relative to the preregistered hypotheses."),
    ("f9b2d82c60e41e07256725a3ad379cac8616bfbe", "A15", "post hoc plan", "Plan for the microarray and missing-gene extensions of the ablation, with its reading rule, committed before those scores. Post hoc relative to the preregistered hypotheses."),
]

TABLES = [
    ("pog_primary", "results/stage4/confirm/primary.tsv", "POG570 preregistered contrasts. Differences are SA-Z minus BASE-Z."),
    ("aux_primary", "results/stage7/confirm/tables/hypothesis_primary.tsv", "Auxiliary preregistered primary contrasts."),
    ("aux_secondary", "results/stage7/confirm/tables/hypothesis_secondary.tsv", "Auxiliary preregistered secondary contrasts. These were not rejected."),
    ("external_metrics", "results/stage8/external/tables/metrics.tsv", "Pre-specified descriptive metrics. Not a hypothesis test."),
    ("subcohort", "results/stage9/subcohort_metrics.tsv", "Post-hoc subcohort metrics."),
    ("set_metrics", "results/stage9/set_metrics.tsv", "Post-hoc set metrics."),
    ("cup_vs_sa", "results/stage9/cup_correct_sa_wrong.tsv", "Post-hoc counts. Sample identifiers are omitted if present."),
    ("mask_metrics", "results/stage9/mask_metrics.tsv", "Post-hoc platform-mask metrics. Models were not retrained."),
    ("mask_random", "results/stage9/mask_random_summary.tsv", "Post-hoc random-drop summary."),
    ("imvigor_contribution", "results/stage9/imvigor_contribution.tsv", "Post-hoc IMvigor210 contributions. Not a marker panel."),
    ("esca_counts", "results/stage9/esca_histology_counts.tsv", "Post-hoc TCGA esophagus histology counts."),
    ("esca_summary", "results/stage9/esophagus_histology_summary.tsv", "Post-hoc esophagus-center assignment counts."),
    ("leave_one_cohort", "results/stage8/diagnostics/leave_one_cohort.tsv", "Post-hoc leave-one-cohort rows. Not a replacement for AH1-AH3."),
]

CBIO = {
    "blca_iatlas_imvigor210_2017",
    "brca_iatlas_anders_2022",
    "mel_dfci_2019",
    "paad_iatlas_prince_2022",
    "prad_su2c_2019__fpkm_capture",
    "prad_su2c_2019__fpkm_polya",
    "prad_fhcrc",
    "skcm_mskcc_2014",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tables() -> None:
    book = Workbook()
    book.remove(book.active)
    for name, rel, description in TABLES:
        frame = pd.read_csv(ROOT / rel, sep="\t")
        id_cols = [c for c in frame.columns if c.lower() in {"sample_id", "patient_id"}]
        frame = frame.drop(columns=id_cols)
        sheet = book.create_sheet(name[:31])
        sheet.append([description])
        sheet.append(list(frame.columns))
        for row in frame.itertuples(index=False):
            sheet.append(list(row))
    dest = BMC / "Additional_file_2.xlsx"
    book.save(dest)
    print(dest, dest.stat().st_size)


def write_history() -> list[str]:
    rows = []
    for prefix, stage, role, note in HISTORY:
        line = subprocess.check_output(
            ["git", "log", "-1", "--format=%H\t%ci\t%s", prefix],
            cwd=ROOT, text=True,
        ).strip()
        full, when, subject = line.split("\t", 2)
        rows.append((full, when, stage, role, subject, note))
    book = Workbook()
    sheet = book.active
    sheet.title = "history"
    sheet.append(["Internal preregistration history. Version control, not a public registry. Times are the git commit times."])
    sheet.append(["commit", "time", "stage", "role", "subject", "note"])
    for row in rows:
        sheet.append(list(row))
    dest = BMC / "Additional_file_3.xlsx"
    book.save(dest)
    print(dest, dest.stat().st_size)
    return [row[0] for row in rows]


def odbl_text() -> str:
    path = ROOT / "config" / "aux_downloads_A5.tsv"
    for line in path.read_text().splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) > 5 and "ODbL" in parts[5]:
            return parts[5]
    raise SystemExit("ODbL sentence not found")


def copy_release() -> None:
    models = RELEASE / "models"
    models.mkdir(parents=True, exist_ok=True)
    for name in ("BASE_Z.joblib", "SA_Z.joblib"):
        shutil.copy2(ROOT / "results" / "stage2" / "models" / name, models / name)
    pool = RELEASE / "host_pool"
    pool.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "results" / "stage5" / "variants" / "gtex_ref_mean_Z.npy", pool / "gtex_ref_mean_Z.npy")
    shutil.copy2(ROOT / "results" / "stage5" / "variants" / "gtex_ref_mean_tissues.json", pool / "gtex_ref_mean_tissues.json")
    (pool / "CITATION.txt").write_text(
        "The data used for the analyses described in this manuscript were obtained from "
        "the GTEx Portal (https://www.gtexportal.org/home/datasets) and from the UCSC Toil "
        "recompute of GTEx (TcgaTargetGtex_rsem_gene_tpm.gz). Please cite the GTEx Consortium "
        "atlas (doi:10.1126/science.aaz1776) and the GTEx Portal download date.\n"
    )
    site = RELEASE / "site_dictionary"
    site.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "scripts" / "mapping_rules.py", site / "mapping_rules.py")
    pred = RELEASE / "predictions"
    if pred.exists():
        shutil.rmtree(pred)
    pred.mkdir()
    shutil.copy2(ROOT / "results" / "stage5" / "posthoc" / "met500_samples.tsv", pred / "MET500_samples.tsv")
    frame = pd.read_csv(ROOT / "results" / "stage8" / "external" / "tables" / "sample_predictions.tsv", sep="\t")
    public = frame.loc[frame["analysis_cohort"] != "POG570"].copy()
    if (public["analysis_cohort"] == "POG570").any():
        raise SystemExit("POG570 rows remained")
    public.to_csv(pred / "public_cohort_sample_predictions.tsv", sep="\t", index=False)
    notice = odbl_text() + "\n"
    (pred / "ODBL_NOTICE.txt").write_text(notice)
    aux = ROOT / "results" / "stage5" / "aux" / "cohorts"
    for cohort_dir in sorted(p for p in aux.iterdir() if p.is_dir()):
        src = cohort_dir / "samples.tsv"
        if not src.exists():
            continue
        if "pog" in cohort_dir.name.lower():
            raise SystemExit(f"unexpected POG directory {cohort_dir.name}")
        target_dir = pred / "aux" / cohort_dir.name
        target_dir.mkdir(parents=True)
        shutil.copy2(src, target_dir / "samples.tsv")
        if cohort_dir.name in CBIO:
            (target_dir / "ODBL_NOTICE.txt").write_text(notice)
        else:
            (target_dir / "SOURCE_NOTE.txt").write_text(
                "These rows are predictions, not the expression matrix. "
                "Retrieve the expression series from GEO or its recorded source by accession.\n"
            )
    shutil.copy2(ROOT / "manuscript" / "env_pins.txt", RELEASE / "env_pins.txt")
    (RELEASE / "LICENSE").write_text(
        "[AUTHOR TO CHOOSE]\n\n"
        "Candidate: MIT License.\n"
        "The investigator chooses the license before any public release. "
        "This file is a placeholder and is not a grant of license.\n"
    )
    (models / "MODEL_CARD.md").write_text(model_card())
    (RELEASE / "README.md").write_text(readme())


def model_card() -> str:
    return """# Model card: BASE-Z and SA-Z

Name: HostMix-TOO is the locked SA-Z procedure.

## Training data

Primary-tumor expression from the UCSC Toil recomputation of TCGA, with GTEx normal-tissue profiles used as the host pool. TARGET samples were counted and not used for training. This package does not contain the expression matrices.

## Features and classes

Within-sample rank features. The gene list is the 5000 B0 symbols. BASE-Z has no scaler. SA-Z uses a StandardScaler fit on the training matrix. Both are L2-logistic models, C = 0.1, coefficient shape (32, 5000). Colorectal is the COADREAD label. Lung combines LUAD and LUSC.

## Augmentation

BASE-Z is trained on real tumor profiles only. SA-Z adds synthetic mixtures of those profiles with host-pool profiles. The mixture weights, mixture count, and training seed are the values locked in the stage-2 model bundle. They were not refit for the manuscript.

## Performance summary

Preregistered POG570, SA-Z minus BASE-Z: at-risk host rate −0.204 (n = 378); evaluation-set top-1 0.0742 (n = 512); native-truth top-1 −0.0659 (n = 91), not significant. Six auxiliary RNA-seq cohorts: host rate −0.148 (n = 427); top-1 0.171 (n = 729); native-truth top-1 0.257 (n = 315). The three secondary auxiliary contrasts were not met.

Pre-specified description, not a hypothesis test: SA-Z had the lowest host rate of the four compared methods on MET500, POG570, and the auxiliary set. It did not have the highest top-1 on every cohort.

## Known limits

Hosts outside the training pool were a failed secondary test. Missing microarray genes move accuracy, and a random drop of the same count does not fully reproduce the platform mask. Assay shift, including IMvigor210 esophagus calls, is a different failure from host attraction. POG570 sample-level predictions are not in this package.

## Citation

Cite the archive (https://doi.org/10.5281/zenodo.23205336) and the GTEx sentence in `host_pool/CITATION.txt`. The article is under review; its DOI will be added on publication.
"""


def readme() -> str:
    return """# Release skeleton

This directory is the reuse package of the public repository https://github.com/juwonk-yonsei/hostmix-too, archived at Zenodo (https://doi.org/10.5281/zenodo.23205336).

## Reproduce

1. Set `project_root` and `figure_dir` in `config.yaml`. The file has no thread-count key.
2. Run the scripts in `scripts/` in numeric order, from the project root, with the environment pins in `env_pins.txt`.
3. Figure scripts `34_figures_A8.py` and `38_figures_A9.py` write to `figure_dir`. They do not retrain models.
4. Locked confirmation models are `models/BASE_Z.joblib` and `models/SA_Z.joblib`. See `models/MODEL_CARD.md`.
5. Host-pool profiles are `host_pool/`. Read `host_pool/CITATION.txt` before reuse.
6. The site dictionary is `site_dictionary/mapping_rules.py`.
7. Sample-level predictions for public cohorts are under `predictions/`. cBioPortal-derived folders include `ODBL_NOTICE.txt`. POG570 sample-level predictions are not included. Aggregate POG570 metrics are only in the manuscript and in Additional file 2.

## Not included

POG570 expression and POG570 sample-level predictions. TCGA and GTEx expression matrices. Those stay at their sources under the terms quoted in the manuscript.
"""


def manifest() -> None:
    lines = []
    for path in sorted(RELEASE.rglob("*")):
        if not path.is_file():
            continue
        if "smoke_out" in path.parts or path.name == "history.bundle":
            continue
        lines.append(f"{sha256(path)}  {path.relative_to(RELEASE)}")
    (RELEASE / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    print("manifest", len(lines))


def main() -> None:
    write_tables()
    commits = write_history()
    copy_release()
    bundle = RELEASE / "history.bundle"
    if bundle.exists():
        bundle.unlink()
    subprocess.check_call(["git", "bundle", "create", str(bundle), "--all"], cwd=ROOT)
    heads = subprocess.check_output(["git", "bundle", "list-heads", str(bundle)], cwd=ROOT, text=True)
    missing = [c for c in commits if c not in heads and not any(h.startswith(c) for h in heads.split())]
    # list-heads prints full hashes; confirm with cat-file via verify
    verify = subprocess.check_output(["git", "bundle", "verify", str(bundle)], cwd=ROOT, text=True, stderr=subprocess.STDOUT)
    print(verify)
    contained = subprocess.check_output(["git", "rev-list", "--all"], cwd=ROOT, text=True).split()
    absent = [c for c in commits if c not in contained]
    print("history commits absent from repo", absent)
    print("bundle bytes", bundle.stat().st_size)
    manifest()
    with (RELEASE / "SHA256SUMS").open("a") as handle:
        handle.write(f"{sha256(bundle)}  history.bundle\n")
    # hashes of additional files
    for name in ("Additional_file_1.docx", "Additional_file_2.xlsx", "Additional_file_3.xlsx", "A_main.docx"):
        path = BMC / name
        if path.exists():
            print(name, path.stat().st_size, sha256(path))


if __name__ == "__main__":
    main()
