"""Draft label maps. Raw strings are copied from the audited metadata.
Unlisted values stay unmapped. COAD and READ share the learning label COADREAD.
"""
from __future__ import annotations

# phenotype column "primary disease or tissue" -> TCGA project code.
DISEASE_TO_PROJECT = {
    "Breast Invasive Carcinoma": "BRCA",
    "Kidney Clear Cell Carcinoma": "KIRC",
    "Lung Adenocarcinoma": "LUAD",
    "Thyroid Carcinoma": "THCA",
    "Head & Neck Squamous Cell Carcinoma": "HNSC",
    "Lung Squamous Cell Carcinoma": "LUSC",
    "Prostate Adenocarcinoma": "PRAD",
    "Brain Lower Grade Glioma": "LGG",
    "Skin Cutaneous Melanoma": "SKCM",
    "Stomach Adenocarcinoma": "STAD",
    "Ovarian Serous Cystadenocarcinoma": "OV",
    "Bladder Urothelial Carcinoma": "BLCA",
    "Liver Hepatocellular Carcinoma": "LIHC",
    "Colon Adenocarcinoma": "COAD",
    "Kidney Papillary Cell Carcinoma": "KIRP",
    "Cervical & Endocervical Cancer": "CESC",
    "Sarcoma": "SARC",
    "Uterine Corpus Endometrioid Carcinoma": "UCEC",
    "Esophageal Carcinoma": "ESCA",
    "Pheochromocytoma & Paraganglioma": "PCPG",
    "Pancreatic Adenocarcinoma": "PAAD",
    "Acute Myeloid Leukemia": "LAML",
    "Glioblastoma Multiforme": "GBM",
    "Testicular Germ Cell Tumor": "TGCT",
    "Thymoma": "THYM",
    "Rectum Adenocarcinoma": "READ",
    "Kidney Chromophobe": "KICH",
    "Mesothelioma": "MESO",
    "Uveal Melanoma": "UVM",
    "Adrenocortical Cancer": "ACC",
    "Uterine Carcinosarcoma": "UCS",
    "Diffuse Large B-Cell Lymphoma": "DLBC",
    "Cholangiocarcinoma": "CHOL",
}

LEARNING_LABEL = {project: ("COADREAD" if project in {"COAD", "READ"} else project) for project in set(DISEASE_TO_PROJECT.values())}

# Learning label -> organ group. Instruction §4.5.
PROJECT_TO_ORGAN = {
    "BRCA": "Breast",
    "PRAD": "Prostate",
    "COADREAD": "Colorectal",
    "LUAD": "Lung",
    "LUSC": "Lung",
    "LIHC": "Liver",
    "CHOL": "Biliary",
    "PAAD": "Pancreas",
    "STAD": "Stomach",
    "ESCA": "Esophagus",
    "HNSC": "HeadNeck",
    "THCA": "Thyroid",
    "KIRC": "Kidney",
    "KIRP": "Kidney",
    "KICH": "Kidney",
    "BLCA": "Bladder",
    "OV": "Ovary",
    "UCEC": "Uterus",
    "UCS": "Uterus",
    "CESC": "Cervix",
    "SKCM": "Melanoma",
    "UVM": "UvealMelanoma",
    "GBM": "Brain",
    "LGG": "Brain",
    "ACC": "Adrenal",
    "PCPG": "Adrenal",
    "TGCT": "Testis",
    "THYM": "Thymus",
    "SARC": "Sarcoma",
    "MESO": "Mesothelioma",
    "DLBC": "Lymphoid",
    "LAML": "Myeloid",
}

# MET500 column `cohort`. MISC and SECR have no TCGA organ and stay unmapped.
MET500_COHORT_TO_ORGAN = {
    "BRCA": "Breast",
    "PRAD": "Prostate",
    "SARC": "Sarcoma",
    "HNSC": "HeadNeck",
    "CHOL": "Biliary",
    "LUNG": "Lung",
    "BLCA": "Bladder",
    "SKCM": "Melanoma",
    "OV": "Ovary",
    "PAAD": "Pancreas",
    "ESCA": "Esophagus",
    "COLO": "Colorectal",
    "STAD": "Stomach",
    "ACC": "Adrenal",
    "KDNY": "Kidney",
    "HCC": "Liver",
    "GBM": "Brain",
    "THCA": "Thyroid",
    "THYM": "Thymus",
    "TGCT": "Testis",
    "MISC": "unmapped",
    "SECR": "unmapped",
}

# POG570 column ANALYSIS_COHORT. Mixed or non-TCGA cohorts stay unmapped.
# CNS-PNS mixes glioma, meningioma, granulosa cell tumor, MPNST, paraganglioma.
# BCC is basal cell carcinoma, which is not the Melanoma organ group.
# PANC includes pancreatic adenocarcinoma and pancreatic neuroendocrine tumor;
# the cohort code is mapped as one label. See the report.
POG570_COHORT_TO_ORGAN = {
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
    "CNS-PNS": "unmapped",
    "SECR": "unmapped",
    "MISC": "unmapped",
    "BCC": "unmapped",
}

# (source, raw_value) -> standard_site.
# source is met500_biopsy_tissue or pog570_biopsy_site.
# Values that are not anatomically explicit are standard_site "other".
SITE_DICTIONARY = {
    ("met500_biopsy_tissue", "liver"): "liver",
    ("met500_biopsy_tissue", "lymph_node"): "lymph_node",
    ("met500_biopsy_tissue", "soft_tissue"): "soft_tissue",
    ("met500_biopsy_tissue", "lung"): "lung",
    ("met500_biopsy_tissue", "bone_marrow"): "bone_marrow",
    ("met500_biopsy_tissue", "skin"): "skin",
    ("met500_biopsy_tissue", "brain"): "brain",
    ("met500_biopsy_tissue", "oral"): "head_neck",
    ("met500_biopsy_tissue", "bladder"): "bladder",
    ("met500_biopsy_tissue", "breast"): "breast",
    ("met500_biopsy_tissue", "thyroid"): "thyroid",
    ("met500_biopsy_tissue", "adrenal"): "adrenal",
    ("met500_biopsy_tissue", "colon"): "colon_rectum",
    ("met500_biopsy_tissue", "cervix"): "cervix",
    ("met500_biopsy_tissue", "prostate"): "prostate",
    ("met500_biopsy_tissue", "pancreas"): "pancreas",
    ("pog570_biopsy_site", "Liver"): "liver",
    ("pog570_biopsy_site", "Liver "): "liver",
    ("pog570_biopsy_site", "Lung"): "lung",
    ("pog570_biopsy_site", "lung"): "lung",
    ("pog570_biopsy_site", "Bronchus"): "lung",
    ("pog570_biopsy_site", "Trachea"): "lung",
    ("pog570_biopsy_site", "Lymph Node"): "lymph_node",
    ("pog570_biopsy_site", "Brain"): "brain",
    ("pog570_biopsy_site", "Intracranial tumor"): "brain",
    ("pog570_biopsy_site", "Pleura"): "pleura",
    ("pog570_biopsy_site", "Omentum"): "omentum",
    ("pog570_biopsy_site", "Omental Mass"): "omentum",
    ("pog570_biopsy_site", "Peritoneum"): "peritoneum",
    ("pog570_biopsy_site", "Peritoneal Mass"): "peritoneum",
    ("pog570_biopsy_site", "Peritoneal Fluid"): "peritoneum",
    ("pog570_biopsy_site", "Skin"): "skin",
    ("pog570_biopsy_site", "Scalp"): "skin",
    ("pog570_biopsy_site", "Breast"): "breast",
    ("pog570_biopsy_site", "Colon"): "colon_rectum",
    ("pog570_biopsy_site", "Rectum"): "colon_rectum",
    ("pog570_biopsy_site", "Ovary"): "ovary",
    ("pog570_biopsy_site", "Pancreas"): "pancreas",
    ("pog570_biopsy_site", "Pancreatic Mass"): "pancreas",
    ("pog570_biopsy_site", "Stomach"): "stomach",
    ("pog570_biopsy_site", "Esophagus"): "esophagus",
    ("pog570_biopsy_site", "Cervix"): "cervix",
    ("pog570_biopsy_site", "Adrenal"): "adrenal",
    ("pog570_biopsy_site", "Bone Marrow"): "bone_marrow",
    ("pog570_biopsy_site", "Spine"): "bone",
    ("pog570_biopsy_site", "Rib"): "bone",
    ("pog570_biopsy_site", "Sacrum "): "bone",
    ("pog570_biopsy_site", "Sternum"): "bone",
    ("pog570_biopsy_site", "Sternal"): "bone",
    ("pog570_biopsy_site", "Scapula Mass"): "bone",
    ("pog570_biopsy_site", "Oral Cavity"): "head_neck",
    ("pog570_biopsy_site", "Tongue"): "head_neck",
    ("pog570_biopsy_site", "Submandibular Gland"): "head_neck",
    ("pog570_biopsy_site", "Sinus"): "head_neck",
    ("pog570_biopsy_site", "Intranasal"): "head_neck",
    ("pog570_biopsy_site", "Maxilla"): "head_neck",
    ("pog570_biopsy_site", "Masticular Space "): "head_neck",
    ("pog570_biopsy_site", "Chest wall"): "soft_tissue",
    ("pog570_biopsy_site", "Abdominal Wall"): "soft_tissue",
    ("pog570_biopsy_site", "Rectus Muscle"): "soft_tissue",
    ("pog570_biopsy_site", "Psoas Muscle"): "soft_tissue",
    ("pog570_biopsy_site", "Gluteus medius"): "soft_tissue",
    ("pog570_biopsy_site", "Arm"): "soft_tissue",
    ("pog570_biopsy_site", "Leg"): "soft_tissue",
    ("pog570_biopsy_site", "Calf"): "soft_tissue",
    ("pog570_biopsy_site", "Thigh Mass"): "soft_tissue",
    ("pog570_biopsy_site", "Buttock"): "soft_tissue",
    ("pog570_biopsy_site", "Flank Mass"): "soft_tissue",
    ("pog570_biopsy_site", "Upper Back Mass"): "soft_tissue",
    ("pog570_biopsy_site", "Lower Back"): "soft_tissue",
    ("pog570_biopsy_site", "Rectroperitoneum"): "soft_tissue",
}

# standard_site -> (gtex tissue names as in phenotype, native projects, native organs)
# GTEx names are the phenotype strings verified in TcgaTargetGTEX_phenotype.
# bone and bone_marrow share the instruction row "bone_marrow, bone".
# skin and subcutaneous share "skin, subcutaneous".
# peritoneum and omentum share "peritoneum, omentum".
# colon_rectum is the instruction row "colon, rectum".
# head_neck is the instruction row "head_neck, oral".
SITE_NATIVE = {
    "liver": (["Liver"], ["LIHC", "CHOL"], ["Liver", "Biliary"]),
    "lung": (["Lung"], ["LUAD", "LUSC", "MESO"], ["Lung", "Mesothelioma"]),
    "lymph_node": ([], ["DLBC"], ["Lymphoid"]),
    "bone_marrow": ([], ["LAML", "DLBC"], ["Myeloid", "Lymphoid"]),
    "bone": ([], ["LAML", "DLBC"], ["Myeloid", "Lymphoid"]),
    "brain": (["Brain - Cortex"], ["GBM", "LGG"], ["Brain"]),
    "adrenal": (["Adrenal Gland"], ["ACC", "PCPG"], ["Adrenal"]),
    "skin": (["Skin - Not Sun Exposed (Suprapubic)"], ["SKCM"], ["Melanoma"]),
    "subcutaneous": (["Skin - Not Sun Exposed (Suprapubic)"], ["SKCM"], ["Melanoma"]),
    "soft_tissue": (["Adipose - Subcutaneous", "Muscle - Skeletal"], ["SARC"], ["Sarcoma"]),
    "kidney": (["Kidney - Cortex"], ["KIRC", "KIRP", "KICH"], ["Kidney"]),
    "pancreas": (["Pancreas"], ["PAAD"], ["Pancreas"]),
    "peritoneum": (["Adipose - Visceral (Omentum)"], ["MESO"], ["Mesothelioma"]),
    "omentum": (["Adipose - Visceral (Omentum)"], ["MESO"], ["Mesothelioma"]),
    "pleura": ([], ["MESO"], ["Mesothelioma"]),
    "bladder": (["Bladder"], ["BLCA"], ["Bladder"]),
    "thyroid": (["Thyroid"], ["THCA"], ["Thyroid"]),
    "ovary": (["Ovary"], ["OV"], ["Ovary"]),
    "breast": (["Breast - Mammary Tissue"], ["BRCA"], ["Breast"]),
    "stomach": (["Stomach"], ["STAD"], ["Stomach"]),
    "colon_rectum": (["Colon - Transverse"], ["COADREAD"], ["Colorectal"]),
    "prostate": (["Prostate"], ["PRAD"], ["Prostate"]),
    "esophagus": (["Esophagus - Mucosa"], ["ESCA"], ["Esophagus"]),
    "head_neck": (["Minor Salivary Gland"], ["HNSC"], ["HeadNeck"]),
    "cervix": (["Cervix - Ectocervix"], ["CESC"], ["Cervix"]),
    "other": ([], [], []),
}

# Simulation host tissue (GTEx phenotype string) -> native organ groups.
# Liver and Lung follow the site rows that name those GTEx tissues.
# Adipose - Subcutaneous is one of the two soft_tissue references (native SARC -> Sarcoma).
SIM_HOST_NATIVE_ORGANS = {
    "Liver": ["Liver", "Biliary"],
    "Lung": ["Lung", "Mesothelioma"],
    "Adipose - Subcutaneous": ["Sarcoma"],
}


def map_site(source: str, raw) -> str:
    if raw is None or (isinstance(raw, float) and raw != raw):
        return "other"
    key = (source, str(raw))
    return SITE_DICTIONARY.get(key, "other")
