#!/usr/bin/env python3
"""Section 7 checks of A16: replacement text in place, removed sentences gone, abstract length,
standardization terms, placeholders, Figure 4 audit and the release package.

Reads the instruction file for the replacement blocks of its section 2 and the pre-A16 copy in
manuscript/archive/bmc_v5 for the removed sentences. Writes manuscript/checks/checks_A16.tsv.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BMC = ROOT / "manuscript" / "bmc"
INSTRUCTION = ROOT.parent / "A_지시서_16.md"
OLD = ROOT / "manuscript" / "archive" / "bmc_v5" / "A_main.md"
SCALE = re.compile(r"\b(scaled|unscaled|scaling)\b", re.IGNORECASE)

INLINE = [
    ("2.2 keywords", "Data augmentation; Feature standardization;"),
    ("2.4 methods", "the HostMix-TOO mixtures without standardization and with C = 0.1 (hereafter the unstandardized mixture model)"),
    ("2.5c figure call", "(Fig. 4a, b)"),
    ("2.8 availability", "the trained baseline and HostMix-TOO models and the unstandardized mixture model of the post hoc ablation"),
    ("3 legend", "**c** Mean loss of top-1 accuracy in MET500 and POG570 when as many randomly chosen classifier genes were removed as each of four microarray platforms lacks (four platforms, ten draws each)."),
]
REMOVED = [
    ("A15 D1", "reproduced most of HostMix-TOO's accuracy loss on microarray profiles"),
    ("A15 D2", "HostMix-TOO's greater sensitivity to missing genes was reproduced by the models with standardized features"),
    ("A15 D3", "Mixtures without feature standardization may therefore transfer better to such profiles"),
    ("old term", "training with mixtures on unscaled features"),
    ("old 2.3", "which part of the training change is responsible and where the remaining errors go"),
    ("old 2.5d", "In the post hoc extension to the microarray layer"),
    ("old 2.6a", "Two results limit the claim. In POG570"),
    ("old 2.6d", "has a plausible explanation: the mixtures teach the classifier"),
    ("old 2.6d", "rather than the mixtures themselves, contributed to the loss"),
    ("old 2.6g", "Several extensions follow directly"),
    ("old 2.7", ". It did not transfer to microarray profiles"),
    ("old 2.9", "rather than to feature scaling or regularization. We report"),
    ("old Fig. 4 legend", "The number above each ablation bar is the share"),
]


def blocks(text: str) -> list[tuple[str, str]]:
    part = text[text.index("## 2. Part A"):text.index("### 2.10")]
    out, heading, current = [], "", []
    for line in part.splitlines() + [""]:
        if line.startswith("### ") or line.startswith("**("):
            heading = line.strip("# *").split("**")[0][:40]
        if line.startswith(">"):
            body = line[1:].strip()
            if body and not re.fullmatch(r"\*\*[A-Za-z]+\*\*", body):
                current.append(body)
        elif current:
            out.append((heading, " ".join(current)))
            current = []
    return out


def words(abstract: str) -> int:
    return len(abstract.split())


def main() -> None:
    main_text = (BMC / "A_main.md").read_text(encoding="utf-8")
    letter = (BMC / "cover_letter.md").read_text(encoding="utf-8")
    supp1 = (BMC / "Additional_file_1.md").read_text(encoding="utf-8")
    flat = re.sub(r"\s+", " ", main_text)
    rows = []
    for heading, body in blocks(INSTRUCTION.read_text(encoding="utf-8")):
        if heading.startswith("2.1"):
            continue
        target = letter if heading.startswith("2.9") else flat
        for sentence in re.split(r"(?<=[.;:])\s+(?=[A-Z(])", body):
            rows.append({"check": f"replacement {heading}", "item": sentence[:70],
                         "result": "present" if sentence.lstrip("…") in target else "MISSING"})
    abstract = main_text[main_text.index("## Abstract"):main_text.index("## Keywords") if "## Keywords" in main_text else None]
    instruction_abstract = " ".join(line[1:].strip() for line in INSTRUCTION.read_text(encoding="utf-8")
                                    .split("### 2.1")[1].split("- 단어 수")[0].splitlines() if line.startswith(">"))
    instruction_abstract = re.sub(r"\*\*", "", instruction_abstract)
    body_abstract = re.sub(r"\*\*|^## Abstract", "", abstract).strip()
    rows.append({"check": "2.1 abstract verbatim", "item": "abstract",
                 "result": "present" if re.sub(r"\s+", " ", body_abstract) == re.sub(r"\s+", " ", instruction_abstract).strip() else "DIFFERS"})
    rows.append({"check": "abstract words", "item": "including the three headings", "result": words(body_abstract)})
    for name, phrase in INLINE:
        rows.append({"check": name, "item": phrase[:70], "result": "present" if phrase in flat else "MISSING"})
    for name, phrase in REMOVED:
        target = letter if name == "old 2.9" else flat
        rows.append({"check": f"removed {name}", "item": phrase[:70], "result": "gone" if phrase not in target else "STILL PRESENT"})
    for name, text in (("main text and legends", main_text), ("cover letter", letter), ("Additional file 1", supp1)):
        rows.append({"check": "scaled/unscaled/scaling", "item": name, "result": len(SCALE.findall(text))})
        rows.append({"check": "{{ }} placeholders", "item": name, "result": text.count("{{")})
    numbers = pd.read_csv(ROOT / "manuscript" / "checks" / "numbers_A16.tsv", sep="\t")
    rows.append({"check": "2.10 numbers", "item": f"{len(numbers)} values", "result": int((~numbers["match"].astype(bool)).sum())})
    audit = pd.read_csv(ROOT / "manuscript" / "checks" / "figure_audit.tsv", sep="\t")
    fig4 = audit.loc[audit["figure"].eq("Figure4")].iloc[0]
    rows.append({"check": "Figure 4", "item": "overlaps + outside page + outside axes + text on axes",
                 "result": int(fig4["overlaps"] + fig4["outside_page"] + fig4["outside_axes"] + fig4["text_on_axes"])})
    rows.append({"check": "Figure 4", "item": "minimum font in the PDF (pt)", "result": float(fig4["min_font_pt_pdf"])})
    rows.append({"check": "Figure 4", "item": "size (mm)", "result": f"{fig4['width_mm']} x {fig4['height_mm']}"})
    models = sorted(path.name for path in (ROOT / "release" / "models").glob("*posthoc*"))
    rows.append({"check": "release", "item": "post hoc model files", "result": ", ".join(models)})
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "manuscript" / "checks" / "checks_A16.tsv", sep="\t", index=False)
    print(table.to_string(index=False))
    bad = table["result"].astype(str).isin(["MISSING", "STILL PRESENT", "DIFFERS"])
    print("problems:", int(bad.sum()))


if __name__ == "__main__":
    main()
