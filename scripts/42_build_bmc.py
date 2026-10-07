"""Assemble the BMC Bioinformatics manuscript with citation tokens.

Does not change result numbers. Tokens are replaced by scripts/41_vancouver_refs.py.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DRAFT = ROOT / "manuscript" / "A_manuscript_draft.md"
OUT = ROOT / "manuscript" / "bmc" / "A_main.md"


def section(text: str, name: str) -> str:
    match = re.search(rf"^## {re.escape(name)}\n(.*?)(?=^## |\Z)", text, re.S | re.M)
    if not match:
        raise SystemExit(f"missing section {name}")
    return match.group(1).strip()


def main() -> None:
    raw = DRAFT.read_text()
    text = re.sub(r"<!--.*?-->", "", raw, flags=re.S)
    replacements = [
        ("SCOPE (Grewal et al., JAMA Network Open 2019)", "SCOPE {GREWAL}"),
        ("CUP-AI-Dx (Zhao et al., EBioMedicine 2020)", "CUP-AI-Dx {ZHAO}"),
        ("(Moiso et al., Cancer Discovery 2022)", "{MOISO}"),
        ("(Vibert et al., Journal of Molecular Diagnostics 2021)", "{VIBERT}"),
        ("(Liu et al., Cancers 2024)", "{LIU}"),
        (
            "Vincent et al. (Bioinformatics 2014, doi:10.1093/bioinformatics/btu044)",
            "Vincent et al. {VINCENT}",
        ),
        (
            "He et al. (Scientific Reports 2023, doi:10.1038/s41598-023-42465-8)",
            "He et al. {HE}",
        ),
        ("Vincent et al. 2014", "Vincent et al. {VINCENT}"),
        ("He et al. 2023", "He et al. {HE}"),
        (
            "Pleasance et al., Nature Cancer 2020 (doi:10.1038/s43018-020-0050-6)",
            "Pleasance et al. {PLEASANCE}",
        ),
        (
            "Robinson et al., Nature 2017 (doi:10.1038/nature23306)",
            "Robinson et al. {ROBINSON}",
        ),
        (
            "The GTEx Consortium, Science 2020 (doi:10.1126/science.aaz1776)",
            "The GTEx Consortium {GTEX_PAPER}",
        ),
        (
            "Mariathasan et al., Nature 2018 (doi:10.1038/nature25501)",
            "Mariathasan et al. {MARIATHASAN}",
        ),
        (
            "Abida et al., PNAS 2019 (doi:10.1073/pnas.1902651116)",
            "Abida et al. {ABIDA}",
        ),
        (
            "Sanghvi et al. (Science Advances 2024)",
            "Sanghvi et al. {SANGHVI}",
        ),
        ("Nagel et al. (iScience 2026)", "Nagel et al. {NAGEL}"),
        ("(doi:10.1016/0197-2456(86)90046-2)", "{DERSIMONIAN}"),
        (
            "from the UCSC Xena public hub",
            "from the UCSC Xena public hub {XENA_MET500}",
        ),
        (
            "Toil TCGA and GTEx matrices used for training are under `data/processed/toil/`.",
            "Toil TCGA and GTEx matrices used for training were downloaded from the UCSC Toil Xena hub {TOIL} and are under `data/processed/toil/`. Open-access GTEx files are described on the GTEx Portal {GTEX_PORTAL}.",
        ),
        ("including GSE50760.", "including GSE50760 {GSE50760}."),
        (
            "GSE209998 and the microarray layer",
            "GSE209998 {GSE209998} and the microarray layer",
        ),
        (
            "cBioPortal datahub files carry the ODC Open Database License statement stored in that table.",
            "cBioPortal datahub files carry the ODC Open Database License statement stored in that table. The study files cited here are IMvigor210 {IMVIGOR}, Anders breast {ANDERS}, DFCI melanoma {MEL}, Prince pancreas {PAAD}, SU2C polyA {SU2C_POLYA}, and SU2C capture {SU2C_CAPTURE}.",
        ),
        (
            "SU2C confirmation files contain polyA library 266, capture library 62",
            "SU2C confirmation files {SU2C_POLYA} {SU2C_CAPTURE} contain polyA library 266, capture library 62",
        ),
        (
            "Esophagus histology uses GDC `diagnoses.primary_diagnosis`",
            "Esophagus histology uses GDC {GDC} `diagnoses.primary_diagnosis`",
        ),
        (
            "The public SCOPE repository used here",
            "The public SCOPE repository {SCOPE_REPO} used here",
        ),
        (
            "CUP-AI-Dx uses `models/inception_net_1d.h5`",
            "CUP-AI-Dx {CUP_REPO} uses `models/inception_net_1d.h5`",
        ),
        ("For SA-Z on GSE41258,", "For SA-Z on GSE41258 {GSE41258},"),
        (
            "genes present on GPL96, GPL20769, GPL15659, or the FHCRC Agilent platform",
            "genes present on GPL96 {GPL96}, GPL20769 {GPL20769}, GPL15659 {GPL15659}, or the FHCRC Agilent platform {FHCRC}",
        ),
    ]
    for old, new in replacements:
        if old not in text:
            raise SystemExit(f"replacement anchor missing: {old[:80]}")
        text = text.replace(old, new, 1)

    background = section(text, "Introduction")
    results = section(text, "Results")
    discussion = section(text, "Discussion")
    methods = section(text, "Methods")
    marker = "Readers who want a single operating rule"
    if marker not in discussion:
        raise SystemExit("conclusion paragraph missing")
    head, conclusions = discussion.split(marker, 1)
    discussion = head.strip()
    conclusions = (marker + conclusions).strip()

    methods = (
        "Paths and extra stored numbers are in Additional file 1. "
        "Stored confirmatory tables are in Additional file 2. "
        "The commit history of the internal preregistration is Additional file 3.\n\n"
        "Internal preregistration, as used here, means that the analysis plan was committed "
        "in version control, with a recorded time, before the corresponding confirmation was "
        "unlocked. It is not registration on a public registry such as the Open Science Framework. "
        "The full commit history is included as Additional file 3.\n\n"
        + methods
        + "\n\n### Use of language-model tools\n\n"
        "[AUTHOR TO CONFIRM] Large language model tools were used to draft the analysis plan "
        "and the stage-wise work instructions (Claude, Anthropic; model name and version "
        "[AUTHOR TO FILL]), and to write and run code, draft the manuscript, and check "
        "references and numbers (Cursor agent; model name and version [AUTHOR TO FILL]). "
        "The confirmatory numbers come from one preregistered run after an internal "
        "preregistration commit. Independent scripts checked Crossref records against the "
        "reference list and checked manuscript numbers against result files. Author "
        "responsibilities for deciding the questions, approving the plan, and reviewing the "
        "manuscript: [AUTHOR TO FILL]. The authors are responsible for the content of this "
        "manuscript. The tools are not authors and are not listed under Authors' contributions "
        "or Acknowledgements."
    )

    abstract = """Background. Classifiers that assign a tissue of origin are usually trained on primary tumors. A metastatic biopsy also contains host tissue, and the tumor transcriptome can move toward the host. A published contamination model mixed liver profiles with liver microRNA and reduced a contamination-related test error. We ask whether a related mixture step, applied to whole-transcriptome ranks and to several host tissues, changes host attraction and accuracy under preregistered tests.

Results. Host attraction is a call of a site-native organ that is not the true organ. On the preregistered POG570 family, the mixture-augmented score minus the unmixed baseline changed the at-risk host rate by −0.204 (n = 378; 95% CI −0.243 to −0.164) and evaluation-set top-1 accuracy by 0.0742 (n = 512; 95% CI 0.0391 to 0.107). The native-truth contrast was −0.0659 (n = 91; 95% CI −0.132 to 0) and was not significant. On six auxiliary RNA-seq cohorts the three preregistered primary contrasts were met, with host-rate difference −0.148 (n = 427), evaluation-set top-1 difference 0.171 (n = 729), and native-truth top-1 difference 0.257 (n = 315). The three secondary contrasts were not met. In a pre-specified descriptive comparison, public classifiers showed host rates from 0.122 to 0.304, and the mixture-augmented score had the lowest host rate on all three analysis cohorts. Top-1 rank was cohort-dependent.

Conclusions. Report top-1 and host rate together, on a written denominator. Mixture augmentation lowered host attraction on the preregistered at-risk sets. That result does not imply that native-truth accuracy rises on every cohort."""

    availability = """Expression matrices used for training and evaluation remain at their sources. Code, the locked BASE-Z and SA-Z model files, the GTEx-derived host-pool profile, the site dictionary, and sample-level predictions for the public cohorts are prepared for release. The repository URL is [TO BE ASSIGNED]. The archive DOI is [TO BE ASSIGNED]. This draft has not been pushed and has not been deposited.

MET500 gene expression, in FPKM, was downloaded from the UCSC Xena public hub {XENA_MET500}. The study is described by Robinson et al. {ROBINSON}. TCGA and GTEx Toil matrices were downloaded from the UCSC Toil Xena hub {TOIL}. The GTEx Consortium atlas is {GTEX_PAPER}. Open-access GTEx files are on the GTEx Portal {GTEX_PORTAL}. The data used for the analyses described in this manuscript were obtained from the GTEx Portal and from the Toil recomputation of GTEx.

POG570 expression was obtained from the BC Cancer Genome Sciences Centre under the POG570 README terms: the data are for research purposes only, and they will not be redistributed without express written permission from the POG coordinating center (poginfo@bcgsc.ca). Sample-level predictions derived from POG570 are not redistributed. Aggregate metrics are reported in this article. The cohort is described by Pleasance et al. {PLEASANCE}.

GEO series were retrieved from GEO, including GSE50760 {GSE50760}, GSE209998 {GSE209998}, and GSE41258 {GSE41258}. Platform files used for the missing-gene description are GPL96 {GPL96}, GPL20769 as read from GSE71729 {GPL20769}, and GPL15659 {GPL15659}.

cBioPortal datahub files used here carry the ODC Open Database License: share and adapt with attribution, keep the derived data open, and share under the same license. The files are IMvigor210 {IMVIGOR}, Anders breast {ANDERS}, DFCI melanoma {MEL}, Prince pancreas {PAAD}, SU2C polyA {SU2C_POLYA}, SU2C capture {SU2C_CAPTURE}, and the FHCRC Agilent matrix {FHCRC}. Sample-level predictions from these files are prepared with that notice.

SCOPE code is the public repository {SCOPE_REPO}. CUP-AI-Dx code is the public repository {CUP_REPO}."""

    legends = """## Figure titles and legends

Figure 1. Stored simulation of host-pull rate.

Each point is the unweighted mean of the stored host-pull rate across tissue-site rows at that mixing weight. In-pool sites are those whose native organs sit in the training host pool. Pool-out sites are the twelve sites named in the text. The panel is a description of the stored simulation table, not a preregistered patient-level test. Source: the stored stage-5 simulation table, summarized in Additional file 2.

Figure 2. Preregistered SA-Z minus BASE-Z differences.

Points are the stored differences. Intervals are the stored patient-bootstrap percentile intervals. POG570 contrasts are H1, H2, and H3. Auxiliary contrasts are AH1–AH3 and the secondary contrasts AS1–AS3. AS contrasts are the stored secondary contrasts. Source: Additional file 2.

Figure 3. Cohort meta-analysis of the primary contrasts.

The vertical line is the stored DerSimonian–Laird estimate. Cohort intervals are the stored effect plus or minus 1.96 times the square root of the stored variance. The seven-cohort RNA-seq set includes GSE209998, which was not part of the six-cohort preregistered family. Source: Additional file 2.

Figure 4. Locked development comparison on MET500.

The heatmap is the stored comparison of representations and model families on the development cohort. It is context for why SA-Z was locked. It is not a second confirmation. Blank cells are stored as missing. Source: Additional file 1.

Figure 5. Public classifiers and post-hoc set splits.

The top row is the pre-specified descriptive comparison on the common-label set. It is not a hypothesis test. Top-1 denominators are 437, 512, and 729. Host-rate denominators are 361, 378, and 427. The bottom row is post hoc and shows native-truth, at-risk, and pool-out at-risk sets. Host rate is not structurally zero on native-truth when a site lists more than one native organ. Source: Additional file 2.

Figure 6. Post-hoc error shift and missing genes.

Panel A counts patients where BASE-Z was wrong and SA-Z was right. Panel B is the fraction of non-esophagus truths called esophagus. Panel C is the GSE41258 colon and rectum slice, selected evaluation n = 183. Panel D is the change in top-1 after restricting genes to a microarray platform and rescoring frozen models. Random-drop means are in Additional file 2 and are not drawn as bars. Source: Additional file 2.

Figure S1. IMvigor210 esophagus-call contributions.

Mean logistic contribution on IMvigor210 bladder native-truth samples predicted as esophagus. BASE-Z n = 145. SA-Z n = 65. The native-truth denominator is 194. The gene list is not a validated marker panel, and exon-union length is not used as a mechanism. Source: Additional file 2.

Figure S2. Esophagus calls and ESCA histology centers.

TCGA-train esophagus n = 145, split into adenocarcinoma 74, squamous carcinoma 65, and other 6. Centers use pure adenocarcinoma and pure squamous samples only. Each selected-evaluation sample called esophagus is assigned by the larger Spearman correlation. Microarray ranks and TCGA centers were not built on the same gene universe. Source: Additional file 2.
"""

    additional = """## Additional files

- File name: Additional file 1
- File format: DOCX (.docx)
- Title of data: Supplementary methods and stored numbers
- Description of data: Paths, preregistration notes, power and pipeline checks, and post-hoc numeric supplements. These rows do not add a test.

- File name: Additional file 2
- File format: XLSX (.xlsx). Microsoft Excel or another spreadsheet program that reads XLSX.
- Title of data: Supplementary tables
- Description of data: One sheet for each stored confirmatory or post-hoc table cited in the text. The first row of each sheet states what the sheet contains. POG570 sample identifiers are not included.

- File name: Additional file 3
- File format: XLSX (.xlsx). Microsoft Excel or another spreadsheet program that reads XLSX.
- Title of data: Internal preregistration history
- Description of data: Commit hash, time, role, and recorded deviation for the stage-3 through stage-8 locks, freezes, unlock, and confirmation outputs. This is version-control history, not a public-registry deposit.
"""

    body = f"""# HostMix-TOO: mixture augmentation and host attraction in metastatic tissue-of-origin classification

## Title page

Title: HostMix-TOO: mixture augmentation and host attraction in metastatic tissue-of-origin classification

Authors: [AUTHOR TO FILL]

Affiliations: [AUTHOR TO FILL]

Corresponding author: [AUTHOR TO FILL]

Email: [AUTHOR TO FILL]

## Abstract

{abstract}

## Keywords

host attraction; cancer of unknown primary; metastasis; tissue of origin; mixture augmentation; preregistration; transcriptome

## Background

{background}

## Methods

{methods}

## Results

{results}

## Discussion

{discussion}

## Conclusions

{conclusions}

## List of abbreviations

BASE-Z: unmixed L2-logistic model on rank features; BLCA: bladder cancer; CI: confidence interval; CUP: cancer of unknown primary; ESCA: esophageal carcinoma; FFPE: formalin-fixed paraffin-embedded; FPKM: fragments per kilobase per million; GDC: Genomic Data Commons; GEO: Gene Expression Omnibus; GTEx: Genotype-Tissue Expression; ODbL: Open Database License; RNA-seq: RNA sequencing; SA-Z: mixture-augmented L2-logistic model on rank features; TCGA: The Cancer Genome Atlas; TPM: transcripts per million.

## Declarations

### Ethics approval and consent to participate

Not applicable. This study used publicly available, de-identified expression matrices and did not recruit participants.

### Consent for publication

Not applicable.

### Availability of data and materials

{availability}

### Competing interests

[AUTHOR TO FILL]

### Funding

[AUTHOR TO FILL]

### Authors' contributions

[AUTHOR TO FILL]

### Acknowledgements

[AUTHOR TO FILL]

## References

<!-- references inserted by scripts/41_vancouver_refs.py -->

{legends}
{additional}
"""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(body)
    print("wrote", OUT, "words_abstract", len(abstract.split()), "words_conclusions", len(conclusions.split()))


if __name__ == "__main__":
    main()
