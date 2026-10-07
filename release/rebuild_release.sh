#!/usr/bin/env bash
# Rebuild the public release package. The default path only prints the steps.
# It does not push, upload, or rewrite the package unless --apply is given,
# and --apply is refused in this copy.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
dry=0
url=""
doi=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) dry=1; shift ;;
    --url) url="$2"; shift 2 ;;
    --doi) doi="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ "$dry" -ne 1 || -z "$url" || -z "$doi" ]]; then
  echo "usage: release/rebuild_release.sh --dry-run --url URL --doi DOI" >&2
  exit 2
fi
models=(
  models/BASE_Z.joblib models/BASE_Z.npz models/BASE_Z.json
  models/SA_Z.joblib models/SA_Z.npz models/SA_Z.json
  models/MIX_Z0_posthoc.joblib models/MIX_Z0_posthoc.npz models/MIX_Z0_posthoc.json
  models/G_symbols.txt models/MODEL_CARD.md
)
echo "dry-run only; the package is not rewritten"
for file in "${models[@]}"; do
  if [[ ! -f "$root/release/$file" ]]; then
    echo "missing model file: release/$file" >&2
    exit 1
  fi
  echo "would include release/$file"
done
echo "would set README and MODEL_CARD URL to $url"
echo "would set README and MODEL_CARD DOI to $doi"
echo "the git history is not public: release/history.bundle stays out of the public snapshot"
echo "would rewrite release/SHA256SUMS (scripts/87_public_snapshot_A19.py sums)"
echo "would build the public snapshot (scripts/87_public_snapshot_A19.py build)"
pog=$(find "$root/release" -iname '*pog*' -o -iname '*POG570*' | grep -v SHA256SUMS || true)
if [[ -n "$pog" ]]; then
  echo "POG570-named paths under release/:"
  echo "$pog"
else
  echo "no POG570-named file under release/"
fi
if find "$root/release/predictions" -iname '*sample*' -print | grep -qi pog; then
  echo "refusing: a POG570 sample-level prediction is in the package" >&2
  exit 1
fi
echo "sample-level POG570 predictions are not in release/predictions"
echo "dry-run finished; nothing was written"
