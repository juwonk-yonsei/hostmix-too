"""POG570 evaluation labels. Metadata only. Expression is not read.

Exception rules run before the cohort table, from rule 1 through rule 7.
A row keeps the first rule that applies. Rows that match none use the cohort table
(rule id "basic"). SECR and BCC have no organ in that table and are exclude.

Keyword search is a casefold substring of TUMOUR_TYPE or HISTOLOGICAL_TYPE.
The token "net" is not searched.

"small cell" is not counted when the match is immediately preceded by "non-" or
"non ". "non-small cell" contains the substring "small cell", and that phrase is
the TCGA non-small-cell group. This exception is locked in the stage-3 preregistration.

Histology codes are exact matches after strip and upper. They are not substrings.
"""
from __future__ import annotations

import pandas as pd

from mapping_rules import SITE_NATIVE, map_site

# Report 1 ANALYSIS_COHORT marginal counts. Used only to build the fake-label
# multinomial. They are not recomputed from per-sample organs.
COHORT_MARGINAL_N = {
    "BRCA": 144,
    "COLO": 87,
    "LUNG": 67,
    "SARC": 47,
    "PANC": 42,
    "OV": 28,
    "CNS-PNS": 19,
    "CHOL": 14,
    "SKCM": 13,
    "SECR": 12,
    "STAD": 11,
    "MISC": 11,
    "UCEC": 11,
    "LYMP": 11,
    "ESCA": 10,
    "HNSC": 7,
    "UVM": 6,
    "KDNY": 5,
    "ACC": 4,
    "CERV": 4,
    "THCA": 4,
    "THYM": 4,
    "PRAD": 3,
    "BCC": 3,
    "HCC": 2,
    "BLCA": 1,
}

# Basic cohort table. CNS-PNS and MISC are absent: only exception rules assign them.
# SECR and BCC are exclude.
BASIC_ORGAN = {
    "BRCA": "Breast",
    "COLO": "Colorectal",
    "LUNG": "Lung",
    "SARC": "Sarcoma",
    "PANC": "Pancreas",
    "OV": "Ovary",
    "CHOL": "Biliary",
    "SKCM": "Melanoma",
    "STAD": "Stomach",
    "UCEC": "Uterus",
    "LYMP": "Lymphoid",
    "ESCA": "Esophagus",
    "HNSC": "HeadNeck",
    "UVM": "UvealMelanoma",
    "KDNY": "Kidney",
    "ACC": "Adrenal",
    "CERV": "Cervix",
    "THCA": "Thyroid",
    "THYM": "Thymus",
    "PRAD": "Prostate",
    "HCC": "Liver",
    "BLCA": "Bladder",
    "SECR": "exclude",
    "BCC": "exclude",
}

NEURO_CODES = {"LUNE", "CSCLC", "SCLC", "ALUCA", "LNET", "PANET", "SCCO"}
LUNG_MESO_CODES = {"PEMESO", "PLMESO", "PLEMESO"}
CNS_BRAIN = {"GBM", "AODG", "ODG", "PAST", "GNG"}
CNS_EXCLUDE = {"MNG", "GRCT", "APE", "MPE", "EPM", "SFTCNS", "PNET"}
HNSC_EXCLUDE = {"NPC", "ODGC", "CSCC", "SNA"}
# MET500 T3 minimum 0.032692468741067175 rounded to 4 decimal places.
BETA_CUT = 0.0327

# Stage-1 cohort counts mapped by the basic table only. Cohorts with no single
# organ (CNS-PNS, MISC, SECR, BCC) are not in this multinomial.
FAKE_ORGAN_WEIGHTS = {
    organ: sum(n for cohort, n in COHORT_MARGINAL_N.items() if BASIC_ORGAN.get(cohort) == organ)
    for organ in sorted({v for v in BASIC_ORGAN.values() if v != "exclude"})
}


def _text(value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value).casefold()


def _code(value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value).strip().upper()


def small_cell_hit(text: str) -> bool:
    needle = "small cell"
    start = 0
    while True:
        found = text.find(needle, start)
        if found < 0:
            return False
        prefix = text[:found]
        if prefix.endswith("non-") or prefix.endswith("non "):
            start = found + len(needle)
            continue
        return True


def keyword_exclude(tumour: str, histology: str) -> bool:
    for word in ("neuroendocrine", "carcinoid"):
        if word in tumour or word in histology:
            return True
    return small_cell_hit(tumour) or small_cell_hit(histology)


def assign_organ(cohort, tumour_type, histological_type) -> tuple[str, str]:
    """Return (organ or exclude or NA, rule id)."""
    cohort_s = "" if cohort is None or (isinstance(cohort, float) and cohort != cohort) else str(cohort).strip()
    tumour = _text(tumour_type)
    histology = _text(histological_type)
    code = _code(histological_type)
    if keyword_exclude(tumour, histology) or code in NEURO_CODES:
        return "exclude", "1"
    if cohort_s == "LUNG" and code in LUNG_MESO_CODES:
        return "Mesothelioma", "2"
    if cohort_s == "CNS-PNS":
        if code in CNS_BRAIN:
            return "Brain", "3"
        if code == "MPNST":
            return "Sarcoma", "3"
        if code == "PGNG":
            return "Adrenal", "3"
        if code in CNS_EXCLUDE:
            return "exclude", "3"
        return "NA", "3-unlisted"
    if cohort_s == "MISC":
        if code == "ESCC":
            return "Esophagus", "4"
        if code == "NSGCT":
            return "Testis", "4"
        return "exclude", "4"
    if cohort_s == "HNSC":
        if code == "OPHSC":
            return "HeadNeck", "5"
        if code == "HNMUCM":
            return "Melanoma", "5"
        if code in HNSC_EXCLUDE:
            return "exclude", "5"
    if cohort_s == "SARC":
        if code in {"ATRT", "MYEC"}:
            return "exclude", "6"
        return "Sarcoma", "6"
    if cohort_s == "OV":
        if code == "OCS":
            return "exclude", "7"
        return "Ovary", "7"
    if cohort_s in BASIC_ORGAN:
        return BASIC_ORGAN[cohort_s], "basic"
    return "NA", "unlisted-cohort"


def label_table(s1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rec in s1.itertuples(index=False):
        organ, rule = assign_organ(rec.ANALYSIS_COHORT, rec.TUMOUR_TYPE, rec.HISTOLOGICAL_TYPE)
        rows.append({
            "PATIENT_ID": rec.PATIENT_ID,
            "ANALYSIS_COHORT": rec.ANALYSIS_COHORT,
            "TUMOUR_TYPE": rec.TUMOUR_TYPE,
            "HISTOLOGICAL_TYPE": rec.HISTOLOGICAL_TYPE,
            "organ": organ,
            "rule": rule,
        })
    return pd.DataFrame(rows)


def native_organs(site: str) -> set[str]:
    tissues, _projects, organs = SITE_NATIVE.get(site, ([], [], []))
    del tissues, _projects
    return set(organs)


def attach_site(s1: pd.DataFrame) -> pd.DataFrame:
    """Biopsy site only. Diagnosis columns are not required."""
    out = s1[["PATIENT_ID", "BIOPSY_SITE", "METASTATIC_OR_RECURRENCE", "TUMOUR_CONTENT"]].copy()
    out["standard_site"] = [map_site("pog570_biopsy_site", value) for value in out["BIOPSY_SITE"]]
    out["TUMOUR_CONTENT_num"] = pd.to_numeric(out["TUMOUR_CONTENT"], errors="coerce")
    return out


def set_flags(organ: str, site: str) -> dict:
    native = native_organs(site)
    in_eval = organ not in {"exclude", "NA"}
    in_native = bool(in_eval and organ in native)
    at_risk = bool(in_eval and native and organ not in native)
    return {
        "in_eval": in_eval,
        "in_native_truth": in_native,
        "at_risk": at_risk,
        "n_native_organs": len(native),
    }
