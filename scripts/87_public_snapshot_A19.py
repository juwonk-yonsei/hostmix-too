#!/usr/bin/env python3
"""Public snapshot of the project repository, without its git history (A19, after the release decision).

The project history contains files that list POG570 patient identifiers with tumor labels or per-patient
results, whose terms of use do not permit redistribution, so the public repository is a snapshot of one project
commit instead of the history. Manuscript drafts, submission notes and superseded check outputs (INTERNAL,
OLD_CHECK) are left out as well.

Usage
  87_public_snapshot_A19.py sums
      rewrite release/SHA256SUMS over the tracked files of release/ (without SHA256SUMS itself, the history
      bundle and __pycache__); run it after the last change to release/ and commit the result.
  87_public_snapshot_A19.py build [--out DIR]
      export HEAD with git archive into DIR (default: the folder hostmix-too-public next to the project),
      remove the WITHHELD paths, write WITHHELD.tsv and SNAPSHOT.txt, check the result and run git init and
      git add in DIR. Nothing is committed or pushed; the commit command is printed. If DIR already holds a
      snapshot repository, its files are replaced and its .git is kept, so a later snapshot becomes a new
      commit of the same public repository.

Checks of build: the tracked files of the project are unchanged against HEAD; no text file of the snapshot
holds two or more distinct POG570 patient identifiers as whole tokens; no text file outside scripts/ holds a
placeholder or the words of the laboratory README (MARKER); no file is over 50 MB; data/, results/
and .venv/ are absent; release/SHA256SUMS matches and covers the release files of the snapshot.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT.parent / "hostmix-too-public"
WITHHELD = {
    "config/pog570_eval_labels.tsv": "POG570 patient identifiers with tumour type and histology from the POG570 "
                                     "release; its SHA-256 is also in config/pog570_eval_labels.sha256",
    "config/prereg_A6.md": "internal plan (Korean) that lists POG570 patient identifiers with per-patient results; "
                           "the English translation config/prereg_A6_en.md is included",
    "reports/": "internal work reports (Korean); some list POG570 patient identifiers",
    "release/history.bundle": "git bundle of the full project history, which contains the files above",
}
INTERNAL = [
    "manuscript/archive/", "manuscript/A_manuscript_draft.md", "manuscript/A_summary_ko.md",
    "manuscript/A_supplement.md", "manuscript/journal_options.md", "manuscript/README.md",
    "manuscript/release_prep.md", "manuscript/bmc/A_main_v3_work.md", "manuscript/bmc/drafts/",
    "manuscript/bmc/review_pack/", "manuscript/bmc/README_ko.md", "manuscript/bmc/reviewer_candidates.md",
    "manuscript/bmc/SUBMISSION_CHECKLIST_ko.md", "manuscript/bmc/submission_checklist.md",
    "manuscript/bmc/cover_letter.md", "manuscript/bmc/cover_letter.docx",
    "manuscript/checks/claims_ledger.md", "manuscript/checks/claims_map.tsv", "manuscript/checks/claims_sample40.md",
    "manuscript/checks/draft_corrections.md", "manuscript/checks/numbers_check.json",
    "manuscript/checks/refs_check.json", "manuscript/checks/vancouver_check.json",
]
OLD_CHECK = re.compile(r"^manuscript/checks/[^/]+_A1[678]\.\w+$")
OLD_CHECK_READ = {
    "manuscript/checks/links_A17.tsv", "manuscript/checks/links_A18.tsv", "manuscript/checks/claims_sample_A17.md",
    "manuscript/checks/claims_sample_A17.tsv", "manuscript/checks/corrections_A17.tsv",
    "manuscript/checks/af3_git_check_A17.tsv", "manuscript/checks/figure_files_A18.tsv",
    "manuscript/checks/figure_fonts_A18.tsv", "manuscript/checks/figure_pdf_props_A18.tsv",
}
MARKER = re.compile(r"\[(?:AUTHOR TO|TO BE ASSIGNED)|Laboratory copy|Not a public release")
TEXT = {".tsv", ".csv", ".txt", ".md", ".json", ".yaml", ".yml", ".py", ".sh", ".cff", ".html", ".R"}
MAX_BYTES = 50 * 1024 ** 2
POG_ID = re.compile(r"(?<![\w.])(?:POG)?(\d{5})(?![\w.])")


def git(*args: str, cwd: Path = ROOT, binary: bool = False):
    out = subprocess.run(["git", *args], cwd=cwd, capture_output=True, check=True)
    return out.stdout if binary else out.stdout.decode()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def tracked(prefix: str = "") -> list[str]:
    return [p for p in git("ls-files", "-z", *([prefix] if prefix else [])).split("\0") if p]


def release_files() -> list[str]:
    return sorted(p for p in tracked("release") if p not in {"release/SHA256SUMS", "release/history.bundle"}
                  and "__pycache__" not in p and "/smoke_out/" not in p)


def sums() -> None:
    lines = [f"{sha((ROOT / p).read_bytes())}  {p[len('release/'):]}" for p in release_files()]
    (ROOT / "release" / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    print(f"release/SHA256SUMS: {len(lines)} files")


def withheld_paths() -> list[str]:
    files = tracked()
    out = []
    for item in WITHHELD:
        hits = [p for p in files if p == item or (item.endswith("/") and p.startswith(item))]
        if not hits:
            raise SystemExit(f"withheld path not tracked: {item}")
        out += hits
    return out


def internal_paths() -> list[str]:
    files = tracked()
    out = []
    for item in INTERNAL:
        hits = [p for p in files if p == item or (item.endswith("/") and p.startswith(item))]
        if not hits:
            raise SystemExit(f"internal path not tracked: {item}")
        out += hits
    return out + [p for p in files if OLD_CHECK.match(p) and p not in OLD_CHECK_READ]


def pog_ids() -> set[str]:
    rows = git("show", "HEAD:config/pog570_eval_labels.tsv").splitlines()[1:]
    return {r.split("\t")[0] for r in rows if r.strip()}


def build(out: Path) -> None:
    changed = [l for l in git("status", "--porcelain", "--untracked-files=no").splitlines() if l.strip()]
    if changed:
        raise SystemExit("tracked files differ from HEAD; commit first:\n" + "\n".join(changed))
    head = git("rev-parse", "HEAD").strip()
    if out.resolve() == ROOT or ROOT in out.resolve().parents:
        raise SystemExit("the snapshot folder must be outside the project")
    if out.exists() and any(out.iterdir()) and not (out / ".git").is_dir():
        raise SystemExit(f"{out} exists, is not empty and is not a snapshot repository")
    out.mkdir(parents=True, exist_ok=True)
    for p in out.iterdir():
        if p.name != ".git":
            shutil.rmtree(p) if p.is_dir() else p.unlink()

    tarfile.open(fileobj=io.BytesIO(git("archive", "--format=tar", "HEAD", binary=True))).extractall(out, filter="data")
    rows = []
    for path in withheld_paths():
        data = git("show", f"HEAD:{path}", binary=True)
        reason = next(r for k, r in WITHHELD.items() if path == k or (k.endswith("/") and path.startswith(k)))
        rows.append(f"{path}\t{sha(data)}\t{len(data)}\t{reason}")
        (out / path).unlink()
    internal = internal_paths()
    for path in internal:
        (out / path).unlink(missing_ok=True)
    for item in [*WITHHELD, *INTERNAL]:
        if item.endswith("/") and (out / item).exists():
            shutil.rmtree(out / item)
    (out / "WITHHELD.tsv").write_text("path\tsha256\tbytes\treason\n" + "\n".join(rows) + "\n")
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    (out / "SNAPSHOT.txt").write_text(
        f"Snapshot of project commit {head}, made {stamp} with scripts/87_public_snapshot_A19.py.\n"
        f"The git history of the project is not included; {len(rows)} tracked files that list POG570 patient "
        "identifiers, or contain such files, are left out (WITHHELD.tsv). A copy of the history without their "
        "contents, in which every commit hash can be verified, is available to editors and reviewers on request.\n"
        f"{len(internal)} further tracked files are not included: manuscript drafts, submission notes and "
        "superseded check outputs.\n")

    problems = []
    ids = pog_ids()
    worst = (0, "")
    files = [p for p in out.rglob("*") if p.is_file() and ".git" not in p.relative_to(out).parts]
    for p in files:
        rel = str(p.relative_to(out))
        if p.stat().st_size > MAX_BYTES:
            problems.append(f"over 50 MB: {rel}")
        if p.suffix in TEXT:
            text = p.read_text(errors="ignore")
            n = len(ids & set(POG_ID.findall(text)))
            worst = max(worst, (n, rel))
            if n >= 2:
                problems.append(f"{n} POG570 patient identifiers: {rel}")
            if not rel.startswith(("scripts/", "release/scripts/")) and MARKER.search(text):
                problems.append(f"placeholder or internal marker: {rel}")
    for name in ("data", "results", ".venv"):
        if (out / name).exists():
            problems.append(f"present: {name}/")
    listed = {}
    for line in (out / "release" / "SHA256SUMS").read_text().splitlines():
        digest, _, name = line.partition("  ")
        listed[f"release/{name}"] = digest
    shipped = {str(p.relative_to(out)) for p in files if str(p.relative_to(out)).startswith("release/")
               and p.name != "SHA256SUMS" and "__pycache__" not in p.parts}
    for name, digest in listed.items():
        if name not in shipped:
            problems.append(f"SHA256SUMS lists a file not in the snapshot: {name}")
        elif sha((out / name).read_bytes()) != digest:
            problems.append(f"SHA256SUMS mismatch: {name}")
    for name in sorted(shipped - set(listed)):
        problems.append(f"release file not in SHA256SUMS: {name}")

    if not (out / ".git").is_dir():
        git("init", "-q", "-b", "main", cwd=out)
    git("add", "-A", cwd=out)
    staged = git("status", "--porcelain", cwd=out).splitlines()
    print(f"snapshot of {head} in {out}")
    print(f"files: {len(files)}; withheld: {len(rows)}; staged changes: {len(staged)}")
    print(f"most POG570 identifiers in one text file: {worst[0]} ({worst[1] or '-'})")
    if problems:
        print("PROBLEMS (do not commit or push):")
        print("\n".join(f"  {p}" for p in problems))
        sys.exit(1)
    print("checks passed. Review, then commit with your own name, for example:")
    print(f'  cd {out} && git -c user.name="NAME" -c user.email="EMAIL" commit -m "Snapshot of {head[:7]}"')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("sums")
    b = sub.add_parser("build")
    b.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    sums() if args.cmd == "sums" else build(args.out)


if __name__ == "__main__":
    main()
