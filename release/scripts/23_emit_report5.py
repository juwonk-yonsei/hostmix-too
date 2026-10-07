"""Write reports/A_보고서_5.md from the stage-5 result files."""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

def _release_config():
    cfg = Path(__file__).resolve().parents[1] / "config.yaml"
    out = {}
    for line in cfg.read_text().splitlines():
        if ":" not in line or line.strip().startswith("#"):
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if key in {"project_root", "figure_dir", "disk_mount"}:
            out[key] = value.strip().strip('"').strip("'")
    return out


def project_root():
    return Path(_release_config()["project_root"])


import pandas as pd

ROOT = project_root()
OUT = ROOT / "reports" / "A_보고서_5.md"

POOL_IN = [
    "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland",
    "Brain - Cortex", "Liver", "Lung", "Muscle - Skeletal",
    "Skin - Not Sun Exposed (Suprapubic)", "Spleen", "Whole Blood",
]
POOL_OUT = [
    "Stomach", "Colon - Transverse", "Pancreas", "Kidney - Cortex",
    "Breast - Mammary Tissue", "Ovary", "Prostate", "Thyroid",
    "Esophagus - Mucosa", "Minor Salivary Gland", "Cervix - Ectocervix", "Bladder",
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def cell(value: object) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", "/").replace("\n", " ")


def md(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    table = frame if columns is None else frame.loc[:, columns]
    headers = [str(name) for name in table.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in table.itertuples(index=False):
        lines.append("| " + " | ".join(cell(value) for value in row) + " |")
    return "\n".join(lines)


def num(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def mean_text(series: pd.Series) -> str:
    return repr(float(series.mean()))


def method_order(methods) -> list[str]:
    first = ["BASE-Z", "SA-Z", "SC-Z", "LD-Z", "NC-Z", "M1-Z"]
    rest = sorted(method for method in methods if method not in first)
    return [method for method in first if method in set(methods)] + rest


def pivot_sim(sim: pd.DataFrame, column: str) -> pd.DataFrame:
    methods = method_order(sim["method"].unique())
    rows = []
    for tissue in POOL_IN + POOL_OUT:
        for rho in (0.6, 0.4):
            sub = sim.loc[(sim["tissue"] == tissue) & (sim["rho"] == rho)]
            row = {"tissue": tissue, "rho": rho}
            for method in methods:
                hit = sub.loc[sub["method"] == method]
                if hit.empty or pd.isna(hit[column].iloc[0]):
                    row[method] = ""
                elif column == "native_truth_top1" and int(hit["n_native_truth"].iloc[0]) == 0:
                    row[method] = ""
                else:
                    row[method] = repr(float(hit[column].iloc[0]))
            rows.append(row)
    return pd.DataFrame(rows)


def sentence_with(abstract: str, needle: str) -> str:
    text = re.sub(r"<[^>]+>", " ", abstract)
    text = re.sub(r"\s+", " ", text)
    index = text.lower().find(needle.lower())
    if index < 0:
        return ""
    start = text.rfind(". ", 0, index)
    start = 0 if start < 0 else start + 2
    end = text.find(". ", index)
    end = len(text) if end < 0 else end + 1
    return text[start:end].strip()


def main() -> None:
    parts: list[str] = []
    add = parts.append
    sim = num(ROOT / "results/stage5/sim_ext/full.tsv")
    paired = read_tsv(ROOT / "results/stage5/tcga_test/paired.tsv")
    sa_base = paired.loc[paired["comparison"] == "SA-Z - BASE-Z"].iloc[0]
    rho = sim.loc[sim["rho"] == 0.6]

    def pool_mean(method: str, tissues: list[str]) -> str:
        rows = rho.loc[(rho["method"] == method) & (rho["tissue"].isin(tissues))]
        if len(rows) != len(tissues):
            raise SystemExit(f"{method} rho 0.6 rows {len(rows)} expected {len(tissues)}")
        return mean_text(rows["host_pull_rate"])

    loho_rows = []
    for tissue in POOL_IN:
        name = f"SA-LOHO-{tissue}"
        row = rho.loc[(rho["method"] == name) & (rho["tissue"] == tissue)]
        if len(row) != 1:
            raise SystemExit(f"missing LOHO row {name}")
        loho_rows.append(float(row["host_pull_rate"].iloc[0]))
    loho_median = float(pd.Series(loho_rows).median())
    base_same = rho.loc[(rho["method"] == "BASE-Z") & (rho["tissue"].isin(POOL_IN)), "host_pull_rate"]
    sa_same = rho.loc[(rho["method"] == "SA-Z") & (rho["tissue"].isin(POOL_IN)), "host_pull_rate"]

    bundles = {
        "RNA-seq polyA": ["prad_su2c_2019__fpkm_polya"],
        "RNA-seq capture": ["prad_su2c_2019__fpkm_capture"],
        "RNA-seq unspecified": [
            "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
            "paad_iatlas_prince_2022", "GSE50760", "GSE209998",
        ],
        "microarray": ["GSE14018", "GSE41258", "GSE71729", "GSE74685", "prad_fhcrc"],
    }
    bundle_lines = []
    for bundle, cohorts in bundles.items():
        n_eval = n_risk = n_native = 0
        for cohort in cohorts:
            samples = read_tsv(ROOT / "results/stage5/aux/cohorts" / cohort / "samples.tsv")
            keep = samples["flag_candidate"] != "True" if "flag_candidate" in samples.columns else pd.Series(True, index=samples.index)
            kept = samples.loc[keep]
            n_eval += int((kept["in_eval"] == "True").sum())
            n_risk += int((kept["in_risk"] == "True").sum())
            n_native += int((kept["in_native"] == "True").sum())
        bundle_lines.append(f"{bundle}: 포함 {len(cohorts)}개, 플래그 제외 후 평가 {n_eval}, 위험 {n_risk}, 고유 장기 정답 {n_native}")

    literature = json.loads((ROOT / "results/stage5/aux/literature_abstracts.json").read_text())
    by_doi = {item["doi"]: item for item in literature}
    status_counts = {"초록만": len(literature), "원문": 0, "미확인": 0}

    add("# A 보고서 5 — SA 변형·기전, POG570 사후 진단, 보조 코호트 수집")
    add("")
    add("## 1. 요약 (수치만, 10줄 이내)")
    add(f"- TCGA-test top-1 BASE-Z {read_tsv(ROOT / 'results/stage5/tcga_test/overall.tsv').set_index('method').loc['BASE-Z', 'top1']}, SA-Z {read_tsv(ROOT / 'results/stage5/tcga_test/overall.tsv').set_index('method').loc['SA-Z', 'top1']}, 차이(SA-Z−BASE-Z) {sa_base['diff']} [ {sa_base['ci_low']}, {sa_base['ci_high']} ]")
    add(f"- 확장 시뮬레이션 ρ=0.6, 조직 가중 없는 평균 host_pull_rate. 풀 안 BASE-Z {pool_mean('BASE-Z', POOL_IN)}, SA-Z {pool_mean('SA-Z', POOL_IN)}. 풀 밖 BASE-Z {pool_mean('BASE-Z', POOL_OUT)}, SA-Z {pool_mean('SA-Z', POOL_OUT)}.")
    add(f"- SA-LOHO 10개 모델의 뺀 조직 ρ=0.6 host_pull_rate 중앙값 {repr(loho_median)}. 같은 10개 조직의 BASE-Z 중앙값 {repr(float(base_same.median()))}, SA-Z 중앙값 {repr(float(sa_same.median()))}.")
    add("- 보조 코호트, 플래그 샘플을 뺀 뒤: " + "; ".join(bundle_lines) + ".")
    add(f"- 문헌 행 확인 상태: 원문 {status_counts['원문']}, 초록만 {status_counts['초록만']}, 미확인 {status_counts['미확인']}.")
    add("")
    add("## 2. 환경과 자원")
    resources = read_tsv(ROOT / "logs/resource_usage.tsv")
    add(md(resources.loc[resources["step"].str.startswith(("14_", "15_", "16_", "18_", "19_", "20_", "21_", "22_"))]))
    add("- `17_aux_screen.py`는 `/usr/bin/time` 기록이 없다.")
    add("- 5단계 스크립트 `stage5_lib.py`, `14_stage5_train.py`, `15_stage5_eval.py`, `16_stage5_post.py`, `17_aux_screen.py`, `18_aux_meta.py`, `19_aux_build.py`, `20_aux_geo.py`, `21_aux_twochannel.py`, `22_aux_fhcrc.py`에서 `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `NUMBA_NUM_THREADS`, worker 상한, taskset 문자열을 검색했고 일치가 없다.")
    add("- 공유 파일 `scripts/stage2_common.py`의 import 시점 OMP 블록과 `scripts/ks_fast.py`의 NUMBA 스레드 상한을 지웠다. 커밋 b954cd2.")
    add("- `scripts/12_confirm.py`와 `scripts/13_posthoc.py`는 스레드 16 설정이 남아 있다. 두 파일은 고치지 않았고 다시 실행하지 않았다.")
    add("- `config/paths.yaml`의 `n_threads: 16`은 5단계 스크립트가 읽지 않았다.")
    _mount = _release_config()["disk_mount"]
    disk = subprocess.check_output(["df", "-h", _mount], text=True).strip().splitlines()[-1]
    add(f"- 시작 시점 `{_mount}` 가용 741G (`logs/disk_stage5_start.txt`). 보고서 작성 시점 `df -h {_mount}`: {disk}.")
    add("")
    section_tcga(add)
    section_simulation(add, sim)
    section_variants(add, sim)
    section_mechanism(add)
    section_posthoc(add)
    section_aux(add)
    section_literature(add, by_doi)
    add("## 10. 잠금 준수 확인")
    add("- `scripts/17_aux_screen.py`, `18_aux_meta.py`, `19_aux_build.py`, `20_aux_geo.py`, `21_aux_twochannel.py`, `22_aux_fhcrc.py`에서 `joblib`, `predict_organ`, `beta` 문자열을 검색했고 일치가 없다. 보조 코호트 파일에 모델 확률을 쓰지 않았다.")
    add("- 4단계 판정 규칙. H1 significant=True, diff=-0.2037037037037037. H2 significant=True, diff=0.07421875. H3 onesided_low=-0.1208791208791209, 기준 > -0.10, significant=False. 4단계 판단 사항은 없다. 보고서 4를 그대로 확정한다.")
    add("- 보조 코호트 사전등록, 동결, 가짜 라벨 절차, 원고 그림, GEO 다음 단계 시리즈, 4단계 §6 그림은 이 단계에서 하지 않았다.")
    add("")
    section_judgments(add)
    section_deviations(add)
    add("## 13. 산출물 목록과 마지막 커밋 해시")
    add("- 자격 기준 커밋 08e70dbf1d70dfc93900881b731bf2baaf7577ce.")
    add("- 스레드 상한 제거와 5단계 학습·평가·사후 스크립트 커밋 b954cd2.")
    add("- 선별·메타데이터 스크립트 커밋 7946fd2.")
    add("- 특징 빌더와 평가 인덱스 수정 커밋 56d180c57534f7c3f15123baa95f86e77843ce1a.")
    add("- 사후 진단 장기 인덱스 수정 커밋 9c5a86781d47269e2a5fbe642f44b3c23dd95d4f.")
    add("- 보고서 본문 커밋: PENDING_BODY_HASH. 이 해시를 본문에 적은 커밋이 그 다음 HEAD이다.")
    add("")
    OUT.write_text("\n".join(parts) + "\n")
    print(OUT, "lines", OUT.read_text().count("\n"))


def section_tcga(add) -> None:
    add("## 3. TCGA-test와 TCGA-met (§2.1)")
    add("### 3.1 전체 표와 짝비교")
    add(md(read_tsv(ROOT / "results/stage5/tcga_test/overall.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/tcga_test/paired.tsv")))
    add("")
    add("- 변형 21개와 BASE-Z, SA-Z. 파일 `results/stage5/tcga_test/variant_overall.tsv`의 본문 표.")
    add(md(read_tsv(ROOT / "results/stage5/tcga_test/variant_overall.tsv")))
    add("")
    add("### 3.2 TCGA-test 장기별 (BASE-Z, SA-Z)")
    add(md(read_tsv(ROOT / "results/stage5/tcga_test/by_organ.tsv")))
    add("")
    add("### 3.3 TCGA-met 표본 유형 묶음별")
    add(md(read_tsv(ROOT / "results/stage5/tcga_met/overall.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/tcga_met/paired.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/tcga_met/by_sample_type.tsv")))
    add("")


def section_simulation(add, sim: pd.DataFrame) -> None:
    add("## 4. 확장 시뮬레이션 (§2.2)")
    add("### 4.1 조직별 GTEx-sim 수, 뺀 조직")
    add("- 풀 안 10개: " + ", ".join(POOL_IN) + ".")
    add("- 풀 밖 12개: " + ", ".join(POOL_OUT) + ".")
    mixes = pd.read_parquet(ROOT / "results/stage5/sim_ext/mixes.parquet")
    counts = mixes.groupby("tissue", as_index=False).agg(n_rows=("host", "size"), n_hosts=("host", "nunique"))
    add(md(counts))
    add("- `results/stage5/sim_ext/dropped_sim_tissues.tsv`는 비어 있다. 기록된 제외 조직은 없다.")
    repro = num(ROOT / "results/stage5/sim_ext/stage2_reproduction.tsv")
    add(f"- 2단계 시뮬레이션 재현 행 수 {len(repro)}. top-1 절대 차이 최댓값 {repr(float(repro['abs_diff'].max()))}.")
    add(md(read_tsv(ROOT / "results/stage5/sim_ext/stage2_reproduction.tsv")))
    add("")
    add("### 4.2 ρ = 0.6과 0.4의 host_pull_rate (행: 조직 × ρ, 열: 방법)")
    add("- `substitute`가 True인 행은 SC-Z가 SA-Z와 같고 LD-Z가 BASE-Z와 같다. 그 표시는 전체 표 파일에 있다.")
    add(md(pivot_sim(sim, "host_pull_rate")))
    add("")
    add("### 4.3 같은 형식의 top-1")
    add(md(pivot_sim(sim, "top1")))
    add("")
    add("### 4.4 같은 형식의 고유 장기 정답 top-1 (분모가 0이면 비움)")
    add(md(pivot_sim(sim, "native_truth_top1")))
    add("")
    add("### 4.5 전체 표 (파일)")
    add("- 파일 `results/stage5/sim_ext/full.tsv`.")
    add("")


def section_variants(add, sim: pd.DataFrame) -> None:
    add("## 5. SA 변형 (§3)")
    add("### 5.1 학습 기록")
    reproduce = json.loads((ROOT / "results/stage5/variants/sa_reproduce.json").read_text())
    add(f"- SA 혼합 기록 재현: n_mix {reproduce['n_mix']}, feature_max_abs {reproduce['feature_max_abs']}, field_match {reproduce['field_match']}.")
    frozen = json.loads((ROOT / "results/stage5/posthoc/pog_frozen_check.json").read_text())
    add(f"- POG570 특징을 다시 계산한 BASE-Z·SA-Z 확률과 동결 예측의 최대 절대 차이 {repr(frozen['max_abs_proba'])}, top-1 불일치 {int(frozen['top1_mismatch'])}, n {int(frozen['n'])}.")
    log = read_tsv(ROOT / "results/stage5/variants/train_log.tsv")
    add(md(log))
    add("- 에포크별 보류 손실 파일: `results/stage5/variants/mlp_loss_BASE-MLP.tsv`, `results/stage5/variants/mlp_loss_SA-MLP.tsv`.")
    add("- `partial_fit` 도중 predict_proba가 확률 합 1 경고를 냈다. 적합은 지정한 에포크 절차로 끝났다. sklearn `early_stopping`은 켜지 않았다.")
    add("")
    add("### 5.2 변형 요약 표 (한 행에 모델 하나)")
    tcga = num(ROOT / "results/stage5/tcga_test/variant_overall.tsv").set_index("method")
    met = num(ROOT / "results/stage5/posthoc/met500_overall.tsv").set_index("method")
    pog = num(ROOT / "results/stage5/posthoc/pog_overall.tsv").set_index("method")
    rho = sim.loc[sim["rho"] == 0.6]
    rows = []
    for method in list(tcga.index):
        block = rho.loc[rho["method"] == method]
        inn = block.loc[block["tissue"].isin(POOL_IN)]
        out = block.loc[block["tissue"].isin(POOL_OUT)]
        rows.append({
            "method": method,
            "tcga_top1": repr(float(tcga.loc[method, "top1"])),
            "tcga_macro_f1": repr(float(tcga.loc[method, "macro_f1"])),
            "sim06_pool_in_host": mean_text(inn["host_pull_rate"]),
            "sim06_pool_in_top1": mean_text(inn["top1"]),
            "sim06_pool_out_host": mean_text(out["host_pull_rate"]),
            "sim06_pool_out_top1": mean_text(out["top1"]),
            "met500_top1": repr(float(met.loc[method, "top1"])),
            "met500_host_rate": repr(float(met.loc[method, "host_rate"])),
            "met500_native_top1": repr(float(met.loc[method, "native_truth_top1"])),
            "met500_gi_internal": int(met.loc[method, "n_gi_internal"]),
            "met500_stomach_pred": int(met.loc[method, "n_pred_stomach"]),
            "met500_stomach_true": int(met.loc[method, "n_pred_stomach_true"]),
            "pog_top1": repr(float(pog.loc[method, "top1"])),
            "pog_host_rate": repr(float(pog.loc[method, "host_rate"])),
            "pog_native_top1": repr(float(pog.loc[method, "native_truth_top1"])),
            "pog_gi_internal": int(pog.loc[method, "n_gi_internal"]),
            "pog_stomach_pred": int(pog.loc[method, "n_pred_stomach"]),
            "pog_stomach_true": int(pog.loc[method, "n_pred_stomach_true"]),
        })
    add(md(pd.DataFrame(rows)))
    add("")
    add("### 5.3 SA-LOHO 표")
    met_site = num(ROOT / "results/stage5/posthoc/met500_by_site.tsv")
    pog_site = num(ROOT / "results/stage5/posthoc/pog_by_site.tsv")

    def site_rate(frame: pd.DataFrame, site: str, method: str) -> str:
        hit = frame.loc[(frame["site"] == site) & (frame["method"] == method)]
        if hit.empty:
            return ""
        return repr(float(hit["host_rate"].iloc[0]))

    rows = []
    for tissue in POOL_IN:
        name = f"SA-LOHO-{tissue}"
        for rho_value in (0.8, 0.6, 0.4):
            block = sim.loc[(sim["tissue"] == tissue) & (sim["rho"] == rho_value)]
            site = str(block["site"].iloc[0])
            rows.append({
                "left_out": tissue,
                "site": site,
                "rho": rho_value,
                "SA-LOHO": repr(float(block.loc[block["method"] == name, "host_pull_rate"].iloc[0])),
                "BASE-Z": repr(float(block.loc[block["method"] == "BASE-Z", "host_pull_rate"].iloc[0])),
                "SA-Z": repr(float(block.loc[block["method"] == "SA-Z", "host_pull_rate"].iloc[0])),
                "MET500_SA-LOHO": site_rate(met_site, site, name),
                "MET500_BASE-Z": site_rate(met_site, site, "BASE-Z"),
                "MET500_SA-Z": site_rate(met_site, site, "SA-Z"),
                "POG570_SA-LOHO": site_rate(pog_site, site, name),
                "POG570_BASE-Z": site_rate(pog_site, site, "BASE-Z"),
                "POG570_SA-Z": site_rate(pog_site, site, "SA-Z"),
            })
    add(md(pd.DataFrame(rows)))
    add("")
    add("### 5.4 짝비교 (MET500, POG570 사후): 변형 − SA-Z, 변형 − BASE-Z")
    add("- MET500.")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/met500_paired.tsv")))
    add("")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/pog_paired.tsv")))
    add("")
    add("### 5.5 MET500·POG570 부위별 host_rate (모든 모델)")
    add("- MET500.")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/met500_by_site.tsv")))
    add("")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/pog_by_site.tsv")))
    add("")


def section_mechanism(add) -> None:
    add("## 6. 기전 (§4)")
    add("### 6.1 Δ_native, Δ_true (ρ = 0.6은 본문 표, 나머지 ρ는 파일)")
    delta = read_tsv(ROOT / "results/stage5/mechanism/delta_summary.tsv")
    add(md(delta.loc[delta["rho"] == "0.6"] if "rho" in delta.columns else delta))
    add("- 나머지 ρ: 파일 `results/stage5/mechanism/delta_summary.tsv`, `results/stage5/mechanism/delta_values.tsv`.")
    add("")
    add("### 6.2 숙주 방향 (조직마다 세 모델의 상위 5 학습 라벨과 값, L_h)")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/alignment.tsv")))
    add("")
    add("### 6.3 계수 감소와 숙주 방향")
    add("- 짝 조직 Spearman과 22개 조직 중 순위. 대상은 `LABEL_HOST`에 있는 학습 라벨이다.")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/coef_spearman.tsv")))
    add("")
    add("- 같은 라벨의 계수 차이 상위 20개 유전자.")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/coef_top_genes.tsv")))
    add("")
    add("### 6.4 Stomach 쏠림")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/gi_distance.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/gi_distance_mean.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/stomach_sim.tsv")))
    add("")
    add("- MET500.")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/met500_stomach_summary.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/met500_stomach_by_truth.tsv")))
    add("")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/pog_stomach_summary.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/pog_stomach_by_truth.tsv")))
    add("")
    add("### 6.5 Sarcoma → Brain")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/sarcoma_sim.tsv")))
    add("")
    add("- MET500.")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/met500_sarcoma_pred.tsv")))
    add("")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/mechanism/pog_sarcoma_pred.tsv")))
    add("")


def section_posthoc(add) -> None:
    add("## 7. POG570 사후 진단과 선택 규칙 후보 (§5)")
    frozen = json.loads((ROOT / "results/stage5/posthoc/pog_frozen_check.json").read_text())
    add(f"- POG570 특징을 다시 계산한 BASE-Z·SA-Z 확률과 동결 예측의 최대 절대 차이 {repr(frozen['max_abs_proba'])}, top-1 불일치 {int(frozen['top1_mismatch'])}, n {int(frozen['n'])}.")
    add("")
    add("### 7.1 고유 장기 정답 불일치 사례 (사후)")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/native_discordant.tsv")))
    add("")
    add("### 7.2 오류와 conformal 집합 (MET500 개발, POG570 사후)")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/conformal_errors.tsv")))
    add("")
    add("### 7.3 선택 규칙 후보 (MET500 개발, POG570 사후)")
    add("- MET500.")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/met500_rules.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/met500_rules_paired.tsv")))
    add("")
    add("- POG570 (사후).")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/pog_rules.tsv")))
    add("")
    add(md(read_tsv(ROOT / "results/stage5/posthoc/pog_rules_paired.tsv")))
    add("")


def section_aux(add) -> None:
    add("## 8. 보조 코호트 (§6)")
    add("### 8.1 적격 기준 커밋 해시와 시각")
    add("- 08e70dbf1d70dfc93900881b731bf2baaf7577ce, 2026-10-01 22:49:08 +0900. 파일 `config/aux_eligibility_A5.md`.")
    add("")
    add("### 8.2 선별 요약 (경로별 선별 수, 포함 수, 제외 기준별 수)과 포함 코호트 표")
    log = read_tsv(ROOT / "results/stage5/aux/screening_log.tsv")
    latest = log.groupby("id", as_index=False).tail(1)
    add(f"- 기록 행 {len(log)}. id별 마지막 행 {len(latest)}.")
    counts = latest["decision"].value_counts()
    for decision, count in counts.items():
        add(f"- 마지막 판정 {decision}: {int(count)}")
    excluded = latest.loc[latest["decision"] == "제외", "failed_criterion"].value_counts()
    for reason, count in excluded.items():
        add(f"- 마지막 제외 기준 {reason}: {int(count)}")
    add("- 경로 × 마지막 판정.")
    add(md(pd.crosstab(latest["path"], latest["decision"]).reset_index()))
    add("")
    add("- 특징 파일을 만든 코호트.")
    features = read_tsv(ROOT / "config/aux_features_A5.tsv")
    add(md(features))
    add("")
    add("### 8.3 코호트별 다운로드와 이용 조건 원문")
    downloads = read_tsv(ROOT / "config/aux_downloads_A5.tsv")
    add(md(downloads))
    add("")
    add("### 8.4 정답 장기 규칙 적용 결과와 제외 사유별 수")
    for cohort_dir in sorted((ROOT / "results/stage5/aux/cohorts").iterdir()):
        samples_path = cohort_dir / "samples.tsv"
        if not samples_path.exists():
            continue
        samples = read_tsv(samples_path)
        if "in_eval" not in samples.columns:
            continue
        rules = samples.loc[samples["in_eval"] != "True", "organ_rule"].value_counts()
        rule_text = ", ".join(f"{name}:{int(count)}" for name, count in rules.items())
        add(f"- {cohort_dir.name}: n={len(samples)}, eval={(samples['in_eval']=='True').sum()}, risk={(samples['in_risk']=='True').sum()}, native={(samples['in_native']=='True').sum()}, patients={samples['patient_id'].nunique()}, organ={','.join(sorted(set(samples['organ'])))}, eval 밖 organ_rule [{rule_text}].")
    add("")
    add("### 8.5 새 부위 사전 초안 전체 (원시값 → 표준 부위)")
    add(md(read_tsv(ROOT / "config/mappings/site_dictionary_aux.tsv")))
    add("")
    add("### 8.6 유전자 범위")
    add("- 8.2의 특징 표와 같다. `n_G_intersect`, `n_B0_missing`, `n_ks_sets_zeroed` 열.")
    add("")
    add("### 8.7 중복 검출")
    rna = [
        "blca_iatlas_imvigor210_2017", "brca_iatlas_anders_2022", "mel_dfci_2019",
        "GSE50760", "GSE209998", "prad_su2c_2019__fpkm_polya", "prad_su2c_2019__fpkm_capture",
    ]
    for cohort in rna:
        samples = read_tsv(ROOT / "results/stage5/aux/cohorts" / cohort / "samples.tsv")
        flagged = samples["flag_candidate"] == "True"
        add(f"- {cohort}: 플래그 {int(flagged.sum())}. 위험 집합 플래그 포함 {(samples['in_risk']=='True').sum()}, 플래그 제외 {(samples.loc[~flagged, 'in_risk']=='True').sum()}.")
        add(md(read_tsv(ROOT / "results/stage5/aux/cohorts" / cohort / "duplicate_quantiles.tsv")))
        add("")
        add(md(read_tsv(ROOT / "results/stage5/aux/cohorts" / cohort / "duplicate_pairs.tsv")))
        add("")
    add("- paad_iatlas_prince_2022의 중복 파일은 log2(TPM+1)이 아니다. 지정한 통계의 분위수와 후보 쌍은 없다.")
    add("- skcm_mskcc_2014는 기호 교집합 0으로 제외했다. 그 중복 파일은 지정한 통계가 아니다.")
    add("- 마이크로어레이 코호트는 중복 선별을 하지 않았다.")
    add("")
    add("### 8.8 집합 크기 (코호트 × 부위, 플랫폼 묶음 합계, 환자 수)")
    site_rows = []
    for cohort_dir in sorted((ROOT / "results/stage5/aux/cohorts").iterdir()):
        samples_path = cohort_dir / "samples.tsv"
        if not samples_path.exists():
            continue
        samples = read_tsv(samples_path)
        if "standard_site" not in samples.columns:
            continue
        for site, part in samples.groupby("standard_site", dropna=False):
            site_rows.append({
                "cohort": cohort_dir.name,
                "standard_site": site,
                "n": len(part),
                "n_eval": int((part["in_eval"] == "True").sum()),
                "n_risk": int((part["in_risk"] == "True").sum()),
                "n_native": int((part["in_native"] == "True").sum()),
                "n_patients": int(part["patient_id"].nunique()),
            })
    add(md(pd.DataFrame(site_rows)))
    add("")


def section_literature(add, by_doi: dict) -> None:
    add("## 9. 관련 연구 표 (§7)")
    rows = [
        ("10.1001/jamanetworkopen.2019.2597", "SCOPE", "Grewal", "86% (11%)", "168 were metastatic"),
        ("10.1016/j.jmoldx.2021.07.009", "TransCUPtomics", "TransCUPtomics", "96%", "79%"),
        ("10.1016/j.ebiom.2020.103030", "CUP-AI-Dx", "CUP-AI-Dx", "96.70%", "86.96%"),
        ("10.1007/s40291-023-00650-5", "Tempus", "Michuda", "91%", "91%"),
        ("10.1101/2025.08.08.669066", "TCUP", "Landau", "98.3 %", "86.7 %"),
        ("10.1038/s41598-022-13665-5", "Hong", "Hong J", "97%", ""),
        ("10.1038/s41598-023-42465-8", "He", "He B", "97.50%", "91.09%"),
        ("10.1016/j.jmoldx.2011.04.004", "Erlander", "Erlander", "87%", "78%"),
        ("10.5858/2006-130-465-mcohcu", "Ma", "Ma", "84%", "82%"),
        ("10.1200/jco.2008.17.9762", "Pathwork", "Monzon", "87.8%", "84.5%"),
        ("10.1016/j.tranon.2021.101016", "GPSai", "Abraham", "94%", ""),
        ("10.1038/s41467-022-31666-w", "CUPLR", "Nguyen", "90%", "58%"),
        ("10.1038/s41591-023-02482-6", "OncoNPC", "Moon", "0.942", ""),
        ("10.1126/sciadv.adn0220", "Sanghvi", "Sanghvi", "more transcriptomically similar", "metastases shifted"),
        ("10.1007/s00401-025-02939-7", "Jürgensen", "Jürgensen", "", "9/9"),
        ("10.1038/nature23306", "MET500", "Robinson", "", ""),
        ("10.1038/s43018-020-0050-6", "POG570", "Pleasance", "", ""),
    ]
    table = []
    for doi, method, author_needle, acc_needle, met_needle in rows:
        item = by_doi[doi]
        abstract = item.get("abstract") or ""
        acc = sentence_with(abstract, acc_needle) if acc_needle else ""
        met = sentence_with(abstract, met_needle) if met_needle else ""
        if acc_needle and acc_needle.lower() not in acc.lower():
            acc = ""
        if met_needle and met_needle.lower() not in met.lower():
            met = ""
        table.append({
            "method": method,
            "authors": item.get("authors") or "",
            "year": item.get("year") or "",
            "journal": item.get("journal") or "",
            "doi": doi,
            "train_or_accuracy_sentence": acc,
            "metastatic_or_named_sentence": met,
            "status": "초록만",
        })
    add("- 각 행의 문장은 저장한 초록에서 그대로 뽑았다. 해당 바늘이 초록에 없으면 칸은 비운다.")
    add(md(pd.DataFrame(table)))
    add("")
    extra = json.loads((ROOT / "results/stage5/aux/literature_extra_search.json").read_text())
    too = ("tissue of origin", "cancer of unknown primary")
    mix = ("in silico mixture", "admixture", "normal contamination", "tumor purity", "biopsy site", "metastatic site")
    kept = []
    for hit in extra["hits"]:
        blob = (hit.get("title") or "").lower()
        if any(term in blob for term in too) and any(term in blob for term in mix):
            kept.append(hit.get("doi"))
    add(f"- 추가 검색: {extra['engine']}, {extra['date']}, hitCount {extra['hitCount']}, 저장한 상위 인용 히트 {len(extra['hits'])}개. 제목에 TOO/CUP 용어와 혼합·부위 용어가 모두 있는 히트 {len(kept)}개. 저장 히트에는 초록이 없다. 1127건 전체는 내리지 않았다.")
    add("")


def section_judgments(add) -> None:
    add("## 11. 판단이 필요한 사항 (번호, 각 항목에 선택지)")
    add("1. GSE74685와 prad_fhcrc 표본 번호가 144개 겹친다. 이번 실행은 GSE74685를 그 번호에 두고, cBio에만 있는 27개를 prad_fhcrc로 두었다. 선택지: (1) 이번 실행을 유지 (2) prad_fhcrc 27개를 빼기 (3) GSE74685를 빼고 cBio 유전자 행렬 171개를 쓰기.")
    add("2. SU2C polyA·capture와 GSE74685·prad_fhcrc를 같은 환자 번호로 맞추지 않았다. 선택지: (1) 두 묶음을 따로 두기 (2) 공통 번호가 확인되면 한쪽으로 합치기.")
    add("3. GSE209998에는 환자 칸이 없어 patient_id를 바코드 전체로 두었다. 선택지: (1) 바코드 1개 = 환자 1명 (2) 바코드 접두부를 나중에 환자 번호로 정하기.")
    add("4. GSE71729에는 환자 칸이 없어 patient_id를 GSM으로 두었다. 선택지: (1) GSM 1개 = 환자 1명 (2) 다른 환자 키가 확인되면 바꾸기.")
    add("5. 표준 부위 `other`는 고유 장기 집합이 비어 위험 집합에 넣지 않았다. 선택지: (1) 위험 집합에서 빼기 (2) 위험 집합에 넣기.")
    add("6. PAAD 파일은 log2 upper-quartile 정규화 계수이고 중복 통계는 log2(TPM+1)이 아니다. 선택지: (1) 순위 특징은 유지하고 지정 중복 통계는 미계산으로 두기 (2) 코호트를 빼기.")
    add("7. GSE74685에서 NE group이 비어 있는 75개와 NA 3개는 CHGA 양성으로 보지 않고 남겼다. 선택지: (1) 남기기 (2) 결측을 제외하기.")
    add("8. skcm_mskcc_2014 발현 파일은 Entrez만 있어 기호 교집합이 0이고 위험 집합 9로 제외했다. 선택지: (1) 제외 유지 (2) Entrez를 기호로 바꾼 뒤 다시 세기.")
    add("9. IMvigor210은 cBio TPM 파일을 썼고 IMvigor210CoreBiologies R 패키지는 받지 않았다. 선택지: (1) cBio TPM 유지 (2) R 패키지 공개 파일이 있으면 바꾸기.")
    add("10. 지시서는 CUPPA라고 적었고, DOI 10.1038/s41467-022-31666-w 초록의 분류기 이름은 CUPLR이다. 선택지: (1) 문헌 표의 방법 이름을 CUPLR로 두기 (2) CUPPA라는 이름의 다른 논문을 따로 찾기.")
    add("")


def section_deviations(add) -> None:
    add("## 12. 지시서와 다르게 한 것")
    add("- `12_confirm.py`와 `13_posthoc.py`의 스레드 16 설정은 남아 있다. 다시 실행하지 않았다.")
    add("- `17_aux_screen.py`는 `/usr/bin/time` 없이 실행했다. RSS는 없다.")
    add("- `15_stage5_eval.py`의 성공한 전체 실행은 TCGA-test 변형 top-1을 쓰지 않았다. 같은 함수를 `--tcga-variants-only`로 다시 실행해 `results/stage5/tcga_test/variant_overall.tsv`를 썼다.")
    add("- cBio 연구 blca_iatlas_imvigor210_2017, brca_iatlas_anders_2022, mel_dfci_2019, paad_iatlas_prince_2022, prad_su2c_2019의 로그 마지막 판정은 선별이다. 포함으로 덮어쓴 행은 없다. 특징 파일은 있다.")
    add("- `20_aux_geo.py`의 첫 실행은 GPL96 주석 헤더에서 멈췄고, 성공한 다시 실행이 `logs/time_20_aux_geo.txt`를 덮어썼다.")
    add("- `21_aux_twochannel.py`의 첫 실행은 GSE71729 특징을 쓴 뒤 GSE74685 파서에서 멈췄다. 성공한 다시 실행이 `logs/time_21_twochannel.txt`를 덮어썼고, 그 실행은 이미 있는 GSE71729 행렬을 다시 읽지 않았다.")
    add("- GSE74685 특성 행은 표본 사이에 정렬되어 있지 않았다. 첫 라벨의 위험 집합 51은 환자 번호가 부위로 들어간 값이다. 표본 순서대로 라벨만 다시 썼고 Z/K는 다시 계산하지 않았다.")
    add("- PAAD는 음수가 없는 행렬로 1e6 합 척도를 맞췄다. 중복 선별은 log2(그 값+1)이다.")
    add("- GPL96 시리즈 행렬 값은 로그 척도로 보이고, 지수로 되돌리지 않았다. 순위 Z/K는 표본별 양수 배율에 불변이다.")
    add("- 2채널 표본 채널 값도 1e6 합으로 다시 척도를 맞췄다. 순위는 그대로다.")
    add("- prad_fhcrc 공개 파일은 유전자 수준이라 프로브 최댓값 규칙을 적용하지 않았다.")
    add("- GEO 제목만 있고 표본 설계가 없는 히트는 보류로 남겼다. 여섯 시작 시리즈만 특징을 만들었다.")
    add("- cBio 화면 이름 Primary Tumor Site는 채취 부위가 아니라서 그 연구들을 기준 3으로 제외했다.")
    add("- SA-MLP와 BASE-MLP의 partial_fit 도중 확률 합 경고가 있었다.")
    add("- `15_stage5_eval.py` 첫 실행은 표본 인덱스가 지표 딕셔너리로 덮여 멈췄고, 두 번째는 `tcga_met` 디렉터리가 없어 멈췄다. 시간 파일은 성공한 실행이다.")
    add("- `16_stage5_post.py` 첫 실행은 장기 이름 리스트를 배열 인덱스로 써서 멈췄다. 시간 파일은 고친 뒤의 실행이다.")
    add("- 계수 표는 `LABEL_HOST`의 학습 라벨만 포함한다.")
    add("- 추가 문헌 검색은 인용순 25개만 저장했다.")
    add("- Grewal 2019 초록 저장본은 `treatment-resistant cancers an`에서 끊긴다. 표의 문장은 그 저장본까지다.")
    add("")


if __name__ == "__main__":
    main()
