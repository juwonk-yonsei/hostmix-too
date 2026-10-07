"""Additional file 3 against git (A17 section 4.11).

For every row of the history sheet: the commit exists and is an ancestor of HEAD, the time column equals
the committer time (and the author time) printed by git, the subject column equals the git subject, and
the role agrees with the commit subject by the keyword rule ROLE_WORDS (every word of one alternative
must appear in the subject, case-insensitive).

Output: manuscript/checks/af3_git_check_A17.tsv
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AF3 = ROOT / "manuscript" / "bmc" / "Additional_file_3.xlsx"
OUT = ROOT / "manuscript" / "checks" / "af3_git_check_A17.tsv"
ROLE_WORDS = {
    "preregistration": [["preregist"], ["plan", "before"]],
    "freeze": [["freeze"]],
    "confirmation output": [["hash", "confirmation"]],
    "model freeze": [["freeze", "model"]],
    "code freeze": [["freeze", "code"]],
    "output path": [["output", "store"]],
    "unlock": [["unlock"]],
    "deviation fix": [["empty", "site"]],
    "addendum": [["addendum"]],
    "input freeze": [["input", "hash", "before"]],
    "label map": [["label map", "before"]],
    "tool output hash": [["output hash"]],
    "post hoc plan": [["plan", "before"]],
}


def git(*args: str) -> tuple[int, str]:
    run = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    return run.returncode, run.stdout.strip()


def main() -> None:
    sheet = pd.read_excel(AF3, sheet_name="history", header=1, dtype=str)
    rows = []
    for r in sheet.itertuples(index=False):
        code, _ = git("cat-file", "-e", f"{r.commit}^{{commit}}")
        exists = code == 0
        ancestor = git("merge-base", "--is-ancestor", r.commit, "HEAD")[0] == 0 if exists else False
        committer = git("show", "-s", "--format=%ci", r.commit)[1] if exists else ""
        author = git("show", "-s", "--format=%ai", r.commit)[1] if exists else ""
        subject = git("show", "-s", "--format=%s", r.commit)[1] if exists else ""
        words = ROLE_WORDS.get(r.role)
        role_ok = bool(words) and any(all(w in subject.lower() for w in alt) for alt in words)
        rows.append({"commit": r.commit, "stage": r.stage, "role": r.role, "af3_time": r.time,
                     "git_committer_time": committer, "git_author_time": author, "exists": exists,
                     "ancestor_of_head": ancestor, "time_ok": r.time == committer and r.time == author,
                     "af3_subject": r.subject, "git_subject": subject, "subject_ok": r.subject == subject,
                     "role_rule": " | ".join("+".join(alt) for alt in words) if words else "NA", "role_ok": role_ok})
    out = pd.DataFrame(rows)
    out["ok"] = out[["exists", "ancestor_of_head", "time_ok", "subject_ok", "role_ok"]].all(axis=1)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"{len(out)} rows, mismatches {int((~out['ok']).sum())}")
    print(out.loc[~out["ok"], ["commit", "role", "exists", "time_ok", "subject_ok", "role_ok"]].to_string()
          if (~out["ok"]).any() else "all rows agree with git")


if __name__ == "__main__":
    main()
