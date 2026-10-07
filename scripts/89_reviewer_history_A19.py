#!/usr/bin/env python3
"""Copy of the project git history for editors and reviewers, without POG570 contents (A19).

The project history contains files that list POG570 patient identifiers with tumor labels or per-patient
results. The copy keeps every commit, tree and tag object unchanged, so every commit hash can be verified, and
leaves out the blobs (file contents) of those files. Rewriting the history would change the commit hashes.

Usage
  89_reviewer_history_A19.py [--out DIR]
      DIR (default: the folder hostmix-too-reviewer-history next to the project) receives
        history.git/                     bare repository: all commits, trees and tags of the branches and tags,
                                         and all blobs except those in MISSING_BLOBS.tsv
        MISSING_BLOBS.tsv                git blob id, bytes, SHA-256, paths in the history, reason
        masked/config/prereg_A6.md       the auxiliary preregistration with the POG570 identifiers masked
        README.md                        how to verify the copy
      and the archive hostmix-too-reviewer-history_<commit>.tar.gz next to DIR, with its SHA-256.

A blob is left out when one of its paths is a WITHHELD path of 87_public_snapshot_A19.py, lies under data/ or
results/, or when its text (or the text of a member, for zip-based files) holds two or more distinct POG570
patient identifiers as whole tokens (the threshold of the public snapshot; single hits in the history are numbers
such as file sizes and OOXML attribute values). Other binary blobs (PDF, PNG, NumPy, joblib, bundle) are not
scanned. Checks: the copy has the same commits as the project; git fsck
reports exactly the left-out blobs as missing; no kept text blob holds a POG570 identifier; every commit hash in
Additional file 3 is present. Nothing is committed, pushed or sent.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import re
import shutil
import subprocess
import tarfile
import zipfile
from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT.parent / "hostmix-too-reviewer-history"
LABELS = "config/pog570_eval_labels.tsv"
MASKED = "config/prereg_A6.md"
REFS = ("--branches", "--tags")
MIN_IDS = 2

spec = importlib.util.spec_from_file_location("snapshot", ROOT / "scripts" / "87_public_snapshot_A19.py")
snapshot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snapshot)
WITHHELD, POG_ID = snapshot.WITHHELD, snapshot.POG_ID


def git(*args: str, cwd: Path = ROOT, binary: bool = False, stdin: bytes | None = None, check: bool = True):
    out = subprocess.run(["git", *args], cwd=cwd, input=stdin, capture_output=True, check=check)
    return out if not check else (out.stdout if binary else out.stdout.decode())


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pog_ids() -> set[str]:
    ids = set()
    for commit in git("log", *REFS, "--format=%H", "--", LABELS).split():
        rows = git("show", f"{commit}:{LABELS}").splitlines()[1:]
        ids |= {r.split("\t")[0] for r in rows if r.strip()}
    return ids


def blob_paths() -> dict[str, set[str]]:
    paths: dict[str, set[str]] = {}
    for commit in git("rev-list", *REFS).split():
        for entry in git("ls-tree", "-r", "-z", "--full-tree", commit).split("\0"):
            if not entry:
                continue
            meta, path = entry.split("\t", 1)
            _, kind, oid = meta.split()
            if kind == "blob":
                paths.setdefault(oid, set()).add(path)
    return paths


def texts(data: bytes) -> list[str] | None:
    if data[:4] == b"PK\x03\x04":
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                return [z.read(n).decode("utf-8", "ignore") for n in z.namelist()]
        except zipfile.BadZipFile:
            return None
    if b"\0" in data[:8000]:
        return None
    return [data.decode("utf-8", "ignore")]


def hits(data: bytes, ids: set[str]) -> int | None:
    parts = texts(data)
    return None if parts is None else len(ids & {m for t in parts for m in POG_ID.findall(t)})


def withheld(path: str) -> bool:
    return any(path == k or (k.endswith("/") and path.startswith(k)) for k in WITHHELD)


def reason(paths: set[str], n: int | None) -> str:
    if any(withheld(p) for p in paths):
        return "withheld path (WITHHELD.tsv of the public snapshot)"
    if any(p.startswith(("data/", "results/")) for p in paths):
        return "data/ or results/"
    if n and n >= MIN_IDS:
        return f"{n} POG570 patient identifiers in the content"
    return ""


def af3_commits() -> set[str]:
    wb = load_workbook(ROOT / "manuscript" / "bmc" / "Additional_file_3.xlsx", read_only=True)
    found = set()
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            found |= {c for c in row if isinstance(c, str) and re.fullmatch(r"[0-9a-f]{7,40}", c)}
    return found


def readme(head: str, n_commits: int, n_missing: int, masked_oid: str, n_masked: int) -> str:
    return f"""# Git history of the HostMix-TOO project, without POG570 contents

This is a copy of the project history up to commit {head}, for editors and reviewers. It is confidential to the
review and is not part of the public repository.

`history.git` is a bare git repository with all {n_commits} commits, all trees and all tags of the project, as
they were written. It lacks the contents (blobs) of {n_missing} files, listed in `MISSING_BLOBS.tsv`: files that
list POG570 patient identifiers with their tumor labels or per-patient results, the internal work reports and the
git bundle of the history. Their git blob ids are still in the trees, so every tree and commit hash, including
those of the preregistration commits in Additional file 3, is the original one.

Verify

    git -C history.git fsck --full --no-dangling
        reports each blob of MISSING_BLOBS.tsv as missing (with the trees that point to it) and nothing else; git checks the hash of every object
        it reads, so the commits and trees present are the original ones
    git -C history.git log --format='%H %ci %s'
        the commit list of Additional file 3
    git -C history.git show --stat COMMIT
    git -C history.git show COMMIT:PATH
        works for every file except those in MISSING_BLOBS.tsv

The commits carry the placeholder identity `aa <aa@aa.com>` under which the analysis agent committed (see Use of
large language models in the manuscript); changing it would change every commit hash.

A full checkout is not possible because of the missing blobs; `git show COMMIT:PATH` and `git ls-tree -r COMMIT`
read single files and file lists.

`masked/{MASKED}` is the Korean auxiliary preregistration (git blob {masked_oid}) with its {n_masked} POG570 patient
identifiers replaced by `[POG570 ID]`; nothing else is changed. Its English translation, `config/prereg_A6_en.md`,
is in the history and in the public repository.

Made with `scripts/89_reviewer_history_A19.py`.
"""


def build(out: Path) -> None:
    changed = [l for l in git("status", "--porcelain", "--untracked-files=no").splitlines() if l.strip()]
    if changed:
        raise SystemExit("tracked files differ from HEAD; commit first:\n" + "\n".join(changed))
    if out.resolve() == ROOT or ROOT in out.resolve().parents:
        raise SystemExit("the output folder must be outside the project")
    head = git("rev-parse", "HEAD").strip()
    ids = pog_ids()
    paths = blob_paths()
    objects = [line.split()[0] for line in git("rev-list", "--objects", *REFS).splitlines() if line.strip()]

    missing, unscanned = {}, 0
    for oid, ps in sorted(paths.items()):
        data = git("cat-file", "blob", oid, binary=True)
        n = hits(data, ids)
        unscanned += n is None
        why = reason(ps, n)
        if why:
            missing[oid] = (len(data), sha256(data), ";".join(sorted(ps)), why)
    keep = [oid for oid in objects if oid not in missing]

    if out.exists():
        shutil.rmtree(out)
    repo = out / "history.git"
    repo.mkdir(parents=True)
    git("init", "-q", "--bare", cwd=repo)
    pack = git("pack-objects", "--stdout", "-q", stdin="\n".join(keep).encode() + b"\n", binary=True)
    git("index-pack", "--stdin", stdin=pack, cwd=repo, binary=True)
    refs = [line.split(" ", 1) for line in git("for-each-ref", "--format=%(objectname) %(refname)",
                                                 "refs/heads", "refs/tags").splitlines()]
    for oid, name in refs:
        git("update-ref", name, oid, cwd=repo)
    branch = git("symbolic-ref", "HEAD").strip()
    git("symbolic-ref", "HEAD", branch, cwd=repo)

    (out / "MISSING_BLOBS.tsv").write_text(
        "blob\tbytes\tsha256\tpaths\treason\n"
        + "".join(f"{oid}\t{b}\t{s}\t{p}\t{r}\n" for oid, (b, s, p, r) in sorted(missing.items(), key=lambda x: x[1][2])))
    masked_oid = git("rev-parse", f"HEAD:{MASKED}").strip()
    original = git("show", f"HEAD:{MASKED}")
    masked = POG_ID.sub(lambda m: "[POG570 ID]" if m.group(1) in ids else m.group(0), original)
    target = out / "masked" / MASKED
    target.parent.mkdir(parents=True)
    target.write_text(masked)
    n_masked = sum(1 for m in POG_ID.finditer(original) if m.group(1) in ids)
    n_commits = len(git("rev-list", *REFS).split())
    (out / "README.md").write_text(readme(head, n_commits, len(missing), masked_oid, n_masked))

    problems = []
    copy_commits = set(git("rev-list", *REFS, cwd=repo).split())
    if copy_commits != set(git("rev-list", *REFS).split()):
        problems.append("commit set differs from the project")
    fsck = git("fsck", "--full", "--no-dangling", cwd=repo, check=False)
    reported = set(re.findall(r"missing blob ([0-9a-f]{40})", fsck.stdout.decode() + fsck.stderr.decode()))
    if reported != set(missing):
        problems.append(f"fsck missing blobs {len(reported)} != left out {len(missing)}")
    other = [l for l in (fsck.stdout.decode() + fsck.stderr.decode()).splitlines()
             if l.strip() and "missing blob" not in l and "broken link" not in l and not l.startswith(("to ", "  "))]
    if other:
        problems.append("fsck reports other problems: " + "; ".join(other[:5]))
    leaked = 0
    for oid in keep:
        if git("cat-file", "-t", oid, cwd=repo).strip() == "blob":
            leaked += (hits(git("cat-file", "blob", oid, cwd=repo, binary=True), ids) or 0) >= MIN_IDS
    if leaked:
        problems.append(f"{leaked} kept blobs hold POG570 identifiers")
    af3 = af3_commits()
    absent = [c for c in af3 if git("cat-file", "-e", f"{c}^{{commit}}", cwd=repo, check=False).returncode]
    if absent:
        problems.append(f"Additional file 3 commits absent: {absent[:5]}")
    if POG_ID.search(masked) and ids & set(POG_ID.findall(masked)):
        problems.append("masked prereg_A6.md still holds identifiers")

    archive = out.parent / f"{out.name}_{head[:7]}.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(out, arcname=out.name)
    digest = sha256(archive.read_bytes())
    archive.with_name(archive.name + ".sha256").write_text(f"{digest}  {archive.name}\n")

    print(f"history copy of {head} in {out}")
    print(f"commits {n_commits}; objects kept {len(keep)} of {len(objects)}; blobs left out {len(missing)}; "
          f"binary blobs not scanned {unscanned}; identifiers masked in {MASKED}: {n_masked}")
    print(f"Additional file 3 commits found: {len(af3) - len(absent)} of {len(af3)}")
    print("left out by reason:", {r: sum(v[3] == r for v in missing.values()) for r in {v[3] for v in missing.values()}})
    print(f"archive {archive} ({archive.stat().st_size} bytes) sha256 {digest}")
    if problems:
        print("PROBLEMS (do not send):")
        print("\n".join(f"  {p}" for p in problems))
        raise SystemExit(1)
    print("checks passed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    build(parser.parse_args().out)


if __name__ == "__main__":
    main()
