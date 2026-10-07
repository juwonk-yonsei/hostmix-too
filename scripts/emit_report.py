#!/usr/bin/env python3
"""Write reports/A_보고서_1.md from the stage-1 result files. No new calculations of accuracy."""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def cell(value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    if isinstance(value, float):
        if value.is_integer() and abs(value) < 1e15:
            return str(int(value))
        return f"{value:.6f}"
    return str(value).replace("|", "\\|")


def md(frame: pd.DataFrame) -> str:
    frame = frame.copy()
    header = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "| " + " | ".join("---" for _ in frame.columns) + " |"
    body = ["| " + " | ".join(cell(v) for v in row) + " |" for row in frame.itertuples(index=False, name=None)]
    return "\n".join([header, rule, *body])


def read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t", dtype=str)


def main() -> None:
    summary = json.loads((ROOT / "results/met500_summary.json").read_text())
    cv = json.loads((ROOT / "results/cv_summary.json").read_text())
    internal = json.loads((ROOT / "results/internal_eval.json").read_text())
    ks = json.loads((ROOT / "results/ks_validation.json").read_text())
    genes = json.loads((ROOT / "data/processed/genes/G_summary.json").read_text())
    eval_n = json.loads((ROOT / "results/met500_eval_n.json").read_text())
    pog = json.loads((ROOT / "results/audit/pog570_structure.json").read_text())
    tcga = json.loads((ROOT / "results/audit/tcga_summary.json").read_text())
    met_struct = json.loads((ROOT / "results/audit/met500_structure.json").read_text())
    toil_scale = json.loads((ROOT / "results/audit/toil_scale.json").read_text())
    met_scale = json.loads((ROOT / "results/audit/met500_fpkm_scale.json").read_text())
    readme = (ROOT / "data/raw/POG570/POG570_README.txt").read_text()
    sim = pd.read_csv(ROOT / "results/simulation.tsv", sep="\t")
    liver = sim.loc[(sim.tissue == "Liver") & (sim.rho == 0.4) & (sim.model == "B1")].iloc[0]

    def block(name: str) -> str:
        row = summary[name]
        return (
            f"{name}: n={row['n']}, top-1={row['top1']:.6f} "
            f"(CI {row['top1_ci'][0]:.6f}–{row['top1_ci'][1]:.6f}), "
            f"H={row['n_host']}/{row['n_error']}={row['H']:.6f} "
            f"(CI {row['H_ci'][0]:.6f}–{row['H_ci'][1]:.6f})"
        )

    parts = []
    parts.append(f"""# A 보고서 1 — 데이터 구축, 잠금, 기준 분류기, 숙주 끌림 오류 측정

## 1. 요약 (수치만, 7줄 이내)
- MET500 평가 n={summary['B1']['n']} (sample_source 평가 단위 {eval_n['n_eval_unit']} 중 organ이 unmapped가 아닌 행). B1 장기 top-1={summary['B1']['top1']:.6f} (CI {summary['B1']['top1_ci'][0]:.6f}–{summary['B1']['top1_ci'][1]:.6f}).
- H(B0)={summary['B0']['n_host']}/{summary['B0']['n_error']}={summary['B0']['H']:.6f} (CI {summary['B0']['H_ci'][0]:.6f}–{summary['B0']['H_ci'][1]:.6f}).
- H(B1)={summary['B1']['n_host']}/{summary['B1']['n_error']}={summary['B1']['H']:.6f} (CI {summary['B1']['H_ci'][0]:.6f}–{summary['B1']['H_ci'][1]:.6f}).
- 시뮬레이션 ρ=0.4, Liver 숙주, B1 장기 top-1={liver['top1']:.6f} (n={int(liver['n'])}). ρ는 RNA 분율이다.
- 부트스트랩 2000회, seed 20261001. H가 정의되지 않은 반복은 다섯 조합 모두 0.
- 계획 커밋 704890ea16d523911845ba17d9e43c14e70a9576 (2026-10-01 16:08:26 +0900). 이 커밋은 §10 실행 전이다.

## 2. 환경과 자원
- Python 3.12.4. numpy 1.26.4, pandas 2.2.2, scipy 1.13.1, scikit-learn 1.4.2, pyarrow 14.0.2, h5py 3.11.0, openpyxl 3.1.2, requests 2.32.2, matplotlib 3.8.4, PyYAML 6.0.1, numba 0.59.1.
- CPU: Intel Xeon E5-2697 v4, 72 logical CPUs (2 sockets × 18 cores × 2 threads). 이 분석의 스레드 상한은 16. RAM 상한은 48GB. `free -g` 총량은 251GB.
- 디스크 시작 (`logs/disk_start.txt`, 2026-10-01T14:56:13+09:00): /HDD8T3 7.3T, used 6.0T, avail 883G, 88%.
- 디스크 종료 (`logs/disk_after.txt`): /HDD8T3 7.3T, used 6.1T, avail 830G, 89%.
- 단계별 wall time과 최대 RSS (`/usr/bin/time -v`, Maximum resident set size).

{md(pd.read_csv(ROOT / "logs/resource_usage.tsv", sep="\\t"))}

- 06_measure 최대 RSS 6,910,108 KB. 48GB 상한 아래이다.

## 3. 다운로드
- `.gz`가 HTTP 403인 파일은 비압축 URL을 받았다: `TCGA_GTEX_category.txt`, `gencode.v23.annotation.gene.probemap`, `M.meta.plus.txt`.
- GEO는 series matrix의 `!series_matrix_table_begin` 앞 헤더만 저장했다. 발현표는 저장하지 않았다.
- ONCOfind git commit `bc4a55e63569bd743412530353a35c3609fd17c8`.
- HGNC 파일 다운로드 시각은 download log의 datetime이다 (2026-10-01).
- 크기와 SHA-256은 `logs/download_log.tsv`와 같다.

{md(pd.read_csv(ROOT / "logs/download_log.tsv", sep="\\t")[["source","file","size_bytes","sha256","status"]])}

### POG570 README 전문

```
{readme.rstrip()}
```

### 데이터 이용 정책 페이지 본문
출처: https://bcgsc.ca/data-release-policy-open-access (HTML을 `data/raw/pog570_data_release_policy.html`에 저장. POG570 잠금 디렉터리 밖).

```
Data made available by Canada's Michael Smith Genome Sciences Centre under an open access data release policy grants unrestricted use, distribution and reproduction in any medium, provided the data is properly cited.
Anyone is free to
Copy, distribute, and display the data
Make derivative data sets
Make commercial use of the data
Under the following conditions: Attribution
Canada's Michael Smith Genome Sciences Centre and/or the author(s) must be given credit
For any reuse or distribution, it must be made clear to others what the license terms of this data are.
Any of these conditions can be waived if the authors gives permission.
Unless otherwise specified, the data is released under the Creative Commons v2.5 Attribution license.
```

README의 "Terms and Conditions for Data Use"와 위 페이지 문장은 서로 다르다. 둘 다 원문 그대로 두었다.

## 4. TCGA·GTEx 정리
- phenotype 열: sample, detailed_category, primary disease or tissue, _primary_site, _sample_type, _gender, _study. 인코딩 latin-1 (UTF-8로 읽으면 byte 0xCA에서 실패).
- 행 수 {tcga['n_phenotype']}. _study: TCGA {tcga['n_tcga']}, GTEX {tcga['n_gtex']}, TARGET {tcga['n_target']}. TARGET는 세기만 하고 학습에 쓰지 않았다.
- 학습 후보 {tcga['n_learning']}. 규칙: (바코드 4번째 필드 앞 2자리 01 AND _sample_type == Primary Tumor) OR (03 AND _sample_type == Primary Blood Derived Cancer - Peripheral Blood AND project == LAML).
- code 01이면서 Primary Tumor가 아닌 행 1개: TCGA-23-1023-01, OV, Recurrent Tumor. 학습 후보에서 제외.
- 같은 환자·같은 code의 최대 행 수는 1. 바코드 정렬 중복 제거로 빠진 행은 0.
- 학습 후보 최소 클래스 CHOL={tcga['learning_min_class']}.
- code 06 = {tcga['n_met06']}. code 11 Solid Tissue Normal = {tcga['n_normal11']}. code 11은 세고 학습·평가에 쓰지 않았다.
- 분할 seed 20261001, test_size 0.2, stratify=learning_label. COAD와 READ의 learning_label은 COADREAD. train 7486, test 1872, 합 9358. 환자 집합은 disjoint.
- TCGA-met 392 중 원발이 TCGA-train에 있는 19행은 exclude_from_eval=primary_in_TCGA-train. 평가 n=373. 원발이 TCGA-test에 있는 06은 평가에 남겼다.

학습 후보 프로젝트별 n (train+test):

{md(pd.read_csv(ROOT / "results/audit/tcga_learning_by_project.tsv", sep="\\t"))}

code 06 프로젝트별 n:

{md(pd.read_csv(ROOT / "results/audit/tcga_met06_by_project.tsv", sep="\\t"))}

train/test 학습 라벨별 n:

{md(pd.read_csv(ROOT / "config/split_tcga.tsv", sep="\\t").groupby(["learning_label","split"]).size().unstack(fill_value=0).reset_index())}

GTEx 조직별 샘플·기증자 수. 기증자는 phenotype의 GTEX 샘플에서 하이픈 앞 두 필드. `Cells - Leukemia Cell Line (Cml)` 70개는 ID가 `K-562-...`이고 이 표에서는 기증자 1로 세어져 있다. 분할에서는 `GTEX-`로 시작하지 않는 70개를 기증자로 만들지 않고 `results/audit/gtex_nondonor_ids.tsv`에 제외했다.

{md(pd.read_csv(ROOT / "results/audit/gtex_tissue_counts.tsv", sep="\\t"))}

제외 후 GTEx-ref 3881, GTEx-sim 3911, 합 7792. 조직 안 기증자 겹침 0. 조직마다 기증자를 정렬한 뒤 Generator(20261001)로 섞고, 앞 n//2를 ref, 나머지를 sim. 홀수 1명은 sim. 조직 순서는 정렬 순서이고 분할은 조직 사이에 독립이 아니다 (Generator를 이어서 쓴다).

## 5. MET500
- 열: Sample_id, sample_type, sample_source, dataset, tissue, cohort, run.id, idx, test, tc, biopsy_tissue.
- 행 868. Sample_id 고유 868. run.id == Sample_id 868. sample_type은 tumor 868. dataset은 mctp 868. sample_source 고유 496.
- Sample_id 안 라이브러리 토큰: capt {met_struct['library_capt']}, poly {met_struct['library_poly']}, 그 외 {met_struct['library_other']}.
- sample_source 496 중 capt+poly 둘 다 372, capt만 112, poly만 12. 둘 다 있는 372쌍은 sample_type, dataset, tissue, cohort, test, tc, biopsy_tissue가 모두 일치 (불일치 0).
- 환자 ID 열은 없다. ID 앞부분으로 환자를 만들지 않았다.
- 적용한 평가 단위: sample_source당 1행. poly가 있으면 poly, 없으면 그 라이브러리. 같은 규칙 안에서 Sample_id 정렬 후 첫 행. 결과 496 (poly 384, capt 112).
- test 열: 라이브러리 기준 FALSE 811, TRUE 57. 평가 단위 기준 FALSE 463, TRUE 33. TRUE를 빼지 않았다.
- tc 결측 0. 라이브러리 868행 분위수 0/10/25/50/75/90/100% = 0.25, 0.30, 0.3875, 0.58, 0.78, 0.87, 0.98.
- 평가에 쓴 매핑 437행의 tc 3분위 경계 (qcut, duplicates=drop): {eval_n['tc_tertile_edges']}.
- cohort 빈도 (라이브러리 868):

{md(pd.read_csv(ROOT / "results/audit/met500_cohort.tsv", sep="\\t"))}

- tissue 빈도 (라이브러리 868):

{md(pd.read_csv(ROOT / "results/audit/met500_tissue.tsv", sep="\\t"))}

- biopsy_tissue 빈도 (라이브러리 868):

{md(pd.read_csv(ROOT / "results/audit/met500_biopsy_tissue.tsv", sep="\\t"))}

- 발현 파일 `M.mx.txt.gz`: 유전자가 행, 첫 칸 이름 sample, 유전자 20979, 샘플 868, ID는 버전 있는 Ensembl (예: ENSG00000186092.4). 값은 선형. min {met_scale['raw_min']}, max {met_scale['raw_max']}, 중앙값 {met_scale['raw_median']} ({met_scale['median_method']}). 음수 0, 비유한값 0.
- G 위 FPKM 합으로 TPM = FPKM/sum(FPKM)*1e6. 합이 0이면 0. 그 외 재정규화는 혼합 시뮬레이션에서만 한다.
- `M.mx.log2.txt.gz`는 받았고 TPM 변환에는 쓰지 않았다.

## 6. POG570 (잠금 범위)
- 발현 파일 열 이름 첫 칸 `genes`. 샘플 열 570. 예시 ID 23764. ID 길이 5. 유전자 58051. 예시 유전자 ENSG00000000003 (버전 없음).
- 값: min {pog['raw_min']}, max {pog['raw_max']}, 중앙값 {pog['raw_median']}, 음수 {pog['n_negative_values']}, 비유한값 {pog['n_nonfinite']}, 값 개수 {pog['n_values']}. README 단위 TPM. 발현 열 ID 570개가 Table S1 PATIENT_ID와 같다.
- 발현 값은 구조 요약(min, median, max, 유전자 ID)에서만 읽었고, 진단·생검 부위와 조인하지 않았다. `load_pog570_expression(unlock=False)`는 POG570Locked를 냈다. lock_refused_without_unlock={pog['lock_refused_without_unlock']}.
- Table S1: 570행, 14열. PATIENT_ID, SAMPLE_ID_DNA, SAMPLE_ID_RNA, AGE, GENDER, TUMOUR_TYPE, HISTOLOGICAL_TYPE, BIOPSY_SITE, BIOPSY_COHORT, ANALYSIS_COHORT, PRIMARY_SITE, METASTATIC_OR_RECURRENCE, TUMOUR_CONTENT, EGAD_ID.
- GENDER 원문 빈도 (뒤 공백 유지): `F` 356, `M` 209, `F ` 3, `M ` 2.
- ANALYSIS_COHORT 결측 0. 빈도:

{md(pd.read_csv(ROOT / "results/audit/pog570_s1_ANALYSIS_COHORT.tsv", sep="\\t"))}

- TUMOUR_TYPE 고유 137, HISTOLOGICAL_TYPE 고유 137, BIOPSY_SITE 고유 77, BIOPSY_COHORT 고유 18, PRIMARY_SITE 고유 37이고 NA 82.
- METASTATIC_OR_RECURRENCE: Metastatic 438, Local 122, Unknown 10.
- TUMOUR_CONTENT 결측 0. 분위수 0/10/25/50/75/90/100% = 17, 34, 45, 58, 75, 86, 100.
- BIOPSY_SITE 빈도는 고유값이 50을 넘겨 감사 파일에 상위 30개만 있다. 전체 77×26 교차표는 `results/audit/pog570_biopsy_site_by_analysis_cohort.tsv` (비영 칸 188). 상위 30:

{md(pd.read_csv(ROOT / "results/audit/pog570_s1_BIOPSY_SITE.tsv", sep="\\t"))}

- Table S2: 1987행, 9열, Patient_ID 고유 466. 열: Patient_ID, Drug_name, Treatment_course_number, Class_1, Class_2, Pathway, Therapy_start_date, Therapy_end_or_biopsy_date, On_treatment_at_biopsy.
- On_treatment_at_biopsy: 0이 1598, 1이 389. Class_1 고유 52, Class_2·Pathway·Drug_name은 고유값이 많아 `results/audit/pog570_s2_*.tsv`에 있다. Class_1 파일은 상위 30개에서 잘려 있다.

## 7. 대응표 (초안 전체)
### 7.1 질환→프로젝트

{md(pd.read_csv(ROOT / "config/mappings/tcga_disease_to_project.tsv", sep="\\t"))}

### 7.2 프로젝트→장기

{md(pd.read_csv(ROOT / "config/mappings/project_to_organ.tsv", sep="\\t"))}

### 7.3 MET500 cohort → 장기

{md(pd.read_csv(ROOT / "config/mappings/met500_cohort_to_organ.tsv", sep="\\t"))}

unmapped: MISC, SECR. 평가 단위에서 MISC 39, SECR 20, 합 59. 이 59는 정확도 분모에서 빠졌다.

### 7.4 POG570 ANALYSIS_COHORT → 장기

{md(pd.read_csv(ROOT / "config/mappings/pog570_cohort_to_organ.tsv", sep="\\t"))}

unmapped: CNS-PNS, SECR, MISC, BCC.

### 7.5 생검 부위 사전

{md(pd.read_csv(ROOT / "config/mappings/site_dictionary.tsv", sep="\\t"))}

사전에 없는 POG570 BIOPSY_SITE는 standard_site=other. 개수 합 48. 값과 n: Abdomen 1, Abdominal Mass 12, Diaphragm 1, Epigastric mass 2, Groin 2, Hip 2, Mediastinal Lesion 1, Mesenteric Mass 1, Neck mass 3, Orbital mass 1, Pelvic Mass 6, Pelvis 3, Pericecal Mass 1, Peripheral Blood 1, Perirenal space 1, Skull Base 1, Spleen 1, Sternal mass 2, Supraclavicular Mass 1, Vagina 1, Vaginal Vault 3, `Vaginal Wall ` 1. MET500 biopsy_tissue 중 사전에 없는 값은 0.

### 7.6 표준 부위 → GTEx 조직, 고유 프로젝트, 고유 장기

{md(pd.read_csv(ROOT / "config/mappings/site_native.tsv", sep="\\t"))}

## 8. GEO 후보 메타데이터 감사
발현표는 받지 않았다. series matrix 헤더만 읽었다.

| accession | n | platform | 수집한 필드 |
| --- | --- | --- | --- |
| GSE41258 | 390 | GPL96 | characteristics `tissue` |
| GSE14018 | 36 | GPL96 | characteristics는 키 없이 `Breast cancer` 36. 생검 장기는 `!Sample_source_name_ch1` |
| GSE74685 | 149 | GPL15659 | ch2 characteristics 키 `tumor site` (샘플당 1개) |
| GSE50760 | 54 | GPL11154 | characteristics `tissue`, `ajcc stage` |

GSE41258 `tissue`: Primary Tumor 186, Normal Colon 54, Polyp 48, Liver Metastasis 47, Lung Metastasis 20, Normal Liver 13, Normal Lung 7, Microadenoma 2, Polyp, high grade 1. `anatomic location`은 결측 138을 포함해 Sigmoid colon 61, Ascending colon 47, Cecum 39, Rectosigmoid 38, Descending colon 31, Transverse colon 18, Rectum 15, Splenic Flexure 2, Hepatic Flexure 1. patient id 고유 275.

GSE14018 source_name: Lung 16, Bone 8, Brain 7, Liver 5. characteristics는 전부 `Breast cancer`.

GSE74685 `tumor site` (ch2, 샘플당 하나): LN 69, LUNG 22, LIVER 21, BONE 20, RETROPERITONEAL 4, ADRENAL 4, PERITONEUM 2, KIDNEY 1, APPENDIX 1, SKIN 1, SCROTUM 1, PERITONEAL 1, SPLEEN 1, RENAL 1. ch1 `tissue`는 149개 모두 `reference pool of LNCaP, DU145, PC3, and CWR22 cell lines`. gender는 male 149.

GSE50760 `tissue`: primary colorectal cancer 18, normal-looking surrounding colonic epithelium 18, metastatic colorectal cancer to the liver 18. `ajcc stage`: Stage IV 54.

## 9. 공통 유전자와 KS
- Toil 심볼 {genes['n_toil_symbols']}, MET500 심볼 {genes['n_met500_symbols']}, POG570 심볼 {genes['n_pog570_symbols']}.
- Toil ID {genes['n_toil_ids']}, probemap으로 매핑된 ID {genes['n_toil_ids_mapped']}. MET500 ID {genes['n_met500_ids']}, 매핑 안 된 ID {genes['n_met500_ids_unmapped']}.
- POG570 ID {genes['n_pog570_ids']}, protein-coding 심볼로 연결되지 않은 ID {genes['n_pog570_unmapped_ids']}. 버전을 뗀 ID가 서로 다른 심볼로 가는 경우 {genes['n_unversioned_symbol_conflicts']}.
- |G|={genes['n_G']}. GMT 원본 {genes['n_gmt_raw']}. G 안 구성원 5개 이상인 세트 {genes['n_sets_kept']}.
- 빠른 KS 검증 1000쌍. 최대 절대차 {ks['max_abs_diff']}. 부호 불일치 {ks['n_sign_mismatch']}. 기준(≤1e-6, 부호 불일치 0) 통과. 구현: fast. 0이 많은 MET500 표본 30개의 최소 0 비율 {ks['high_zero_fraction_min']}.
- 저장: `data/processed/features/ks_tcga.parquet` (9750행), `ks_gtex_ref.parquet` (1650행, missing 0), `ks_met500.parquet` (868행), `ks_sim.parquet`. POG570 KS는 없다.

## 10. 기준 분류기
- B0 C={cv['B0']['C']}. 5겹 macro-F1 평균: {cv['B0']['mean_macro_f1']}.
- B1 C={cv['B1']['C']}. 5겹 macro-F1 평균: {cv['B1']['mean_macro_f1']}.
- 동점이면 더 작은 C. 최종 적합 n_iter 최대: B0 109, B1 40. max_iter 5000에 닿지 않았다.
- B0 유전자 필터는 TCGA-train 전체에서 한 번. CV 폴드 안에서 다시 고르지 않았다. B1 StandardScaler는 CV 폴드 안에서 적합하고, 최종 스케일러는 TCGA-train 전체에 적합했다.
- macro-F1의 labels는 그 평가 집합의 y_true에 있는 클래스만이다.

TCGA-test n=1872.

| level | model | top-1 | top-3 | macro-F1 |
| --- | --- | --- | --- | --- |
| project | B0 | {internal['TCGA-test']['project_B0']['top1']:.6f} | {internal['TCGA-test']['project_B0']['top3']:.6f} | {internal['TCGA-test']['project_B0']['macro_f1']:.6f} |
| project | B1 | {internal['TCGA-test']['project_B1']['top1']:.6f} | {internal['TCGA-test']['project_B1']['top3']:.6f} | {internal['TCGA-test']['project_B1']['macro_f1']:.6f} |
| organ | B0 | {internal['TCGA-test']['organ_B0']['top1']:.6f} | {internal['TCGA-test']['organ_B0']['top3']:.6f} | {internal['TCGA-test']['organ_B0']['macro_f1']:.6f} |
| organ | B1 | {internal['TCGA-test']['organ_B1']['top1']:.6f} | {internal['TCGA-test']['organ_B1']['top3']:.6f} | {internal['TCGA-test']['organ_B1']['macro_f1']:.6f} |

TCGA-test 프로젝트별 recall, B0:

{md(pd.DataFrame(internal['TCGA-test']['project_B0']['recall']))}

TCGA-test 프로젝트별 recall, B1:

{md(pd.DataFrame(internal['TCGA-test']['project_B1']['recall']))}

TCGA-test 장기별 recall, B0:

{md(pd.DataFrame(internal['TCGA-test']['organ_B0']['recall']))}

TCGA-test 장기별 recall, B1:

{md(pd.DataFrame(internal['TCGA-test']['organ_B1']['recall']))}

TCGA-met 평가 n=373. train과 환자가 겹쳐 제외된 06은 19.

| level | model | top-1 | top-3 | macro-F1 |
| --- | --- | --- | --- | --- |
| project | B0 | {internal['TCGA-met']['project_B0']['top1']:.6f} | {internal['TCGA-met']['project_B0']['top3']:.6f} | {internal['TCGA-met']['project_B0']['macro_f1']:.6f} |
| project | B1 | {internal['TCGA-met']['project_B1']['top1']:.6f} | {internal['TCGA-met']['project_B1']['top3']:.6f} | {internal['TCGA-met']['project_B1']['macro_f1']:.6f} |
| organ | B0 | {internal['TCGA-met']['organ_B0']['top1']:.6f} | {internal['TCGA-met']['organ_B0']['top3']:.6f} | {internal['TCGA-met']['organ_B0']['macro_f1']:.6f} |
| organ | B1 | {internal['TCGA-met']['organ_B1']['top1']:.6f} | {internal['TCGA-met']['organ_B1']['top3']:.6f} | {internal['TCGA-met']['organ_B1']['macro_f1']:.6f} |

TCGA-met 프로젝트별 top-1:

{md(pd.read_csv(ROOT / "results/tcga_met_by_project.tsv", sep="\\t"))}

SKCM만, Xena `SKCM_clinicalMatrix`의 `tumor_tissue_site`와 15자 sample id를 정확히 조인. 임상 파일 전체 빈도: Regional Lymph Node 226, Primary Tumor 105, Regional Cutaneous or Subcutaneous Tissue (includes satellite and in-transit metastasis) 76, Distant Metastasis 70, 결측 3, [Discrepancy] 1.

{md(pd.read_csv(ROOT / "results/skcm_by_tumor_tissue_site.tsv", sep="\\t"))}

## 11. MET500 결과
평가 n={summary['B0']['n']}. 부트스트랩 단위는 sample_source. 다섯 조합의 n이 같아서 seed 20261001의 복원 추출 인덱스는 같다. NA 예측은 오답이다. M1에서 남은 확률 합이 0이면 예측은 NA.

### 11.1 조합별 정확도

{md(pd.DataFrame([{
    "model": name,
    "n": summary[name]["n"],
    "top1": summary[name]["top1"],
    "top1_ci_low": summary[name]["top1_ci"][0],
    "top1_ci_high": summary[name]["top1_ci"][1],
    "top3": summary[name]["top3"],
    "top3_ci_low": summary[name]["top3_ci"][0],
    "top3_ci_high": summary[name]["top3_ci"][1],
    "macro_f1": summary[name]["macro_f1"],
} for name in ("B0","B0+M1","B1","B1+M1","B1+V0")]))}

### 11.2 오류 분해
H = n_host / n_error. host_rate = n_host / n_at_risk. H CI에서 n_error=0인 반복은 빼고 개수를 센다. 이번 실행에서는 그 개수가 0이다.

{md(pd.DataFrame([{
    "model": name,
    "n_error": summary[name]["n_error"],
    "n_host": summary[name]["n_host"],
    "H": summary[name]["H"],
    "H_ci_low": summary[name]["H_ci"][0],
    "H_ci_high": summary[name]["H_ci"][1],
    "n_at_risk": summary[name]["n_at_risk"],
    "host_rate": summary[name]["host_rate"],
    "host_rate_ci_low": summary[name]["host_rate_ci"][0],
    "host_rate_ci_high": summary[name]["host_rate_ci"][1],
    "n_H_undefined_bootstrap": summary[name]["n_bootstrap_H_undefined"],
} for name in ("B0","B0+M1","B1","B1+M1","B1+V0")]))}

### 11.3 부위·tc 3분위·라이브러리 토큰 (B0, B1)
`factor=library`의 capt/poly는 Sample_id 안 토큰이다. 열 이름 sample_source는 검체 ID이고, 평가 단위로 접은 뒤 행마다 하나씩이다. H가 빈 칸이면 n_error=0이라 H가 정의되지 않는다. host_rate가 빈 칸이면 n_at_risk=0이다.

{md(pd.read_csv(ROOT / "results/met500_error_by_factor.tsv", sep="\\t"))}

### 11.4 생검 부위 고유 장기가 정답인 샘플의 top-1

{md(pd.read_csv(ROOT / "results/met500_native_truth_accuracy.tsv", sep="\\t"))}

### 11.5 B1 혼동 상위 15 (대각 제외)

{md(pd.read_csv(ROOT / "results/b1_confusion_top15.tsv", sep="\\t"))}

전체 칸은 `results/b1_confusion_all.tsv`.

B1 숙주 끌림 오류 51건의 standard_site 수:

{md(pd.read_csv(ROOT / "results/host_errors_B1_by_site.tsv", sep="\\t"))}

다섯 조합의 숙주 끌림 오류 목록(sample id, sample_source, site, tc, truth, pred, top3)은 `results/host_errors_*.tsv`.

### 11.6 장기별 n과 top-1

{md(pd.read_csv(ROOT / "results/met500_by_organ.tsv", sep="\\t"))}

## 12. 시뮬레이션 결과
TCGA-test에서 learning_label 층화로 1000개. seed 20261001. 숙주는 GTEx-sim의 Liver, Lung, Adipose - Subcutaneous. 정규 조직은 복원 추출. 순서: 조직은 위 순서, ρ는 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2. ρ=1에서도 같은 난수 호출로 숙주를 뽑고 (1-ρ)=0을 곱한다.
혼합: G 위에서 각 프로파일을 합 1e6으로 맞춘 뒤 x = ρ·tumor + (1-ρ)·normal. ρ는 RNA 분율이고 세포 분율이 아니다.
host_pull_rate 분모는 정답 장기가 그 조직의 고유 장기 집합 밖인 종양 수. 고유 장기: Liver → Liver, Biliary. Lung → Lung, Mesothelioma. Adipose - Subcutaneous → Sarcoma.
오답 상위 3개 중 빈 칸은 그 순위까지 오답 장기가 없다는 뜻이다.

{md(sim)}

## 13. POG570 위험 대상 수 (부위별)
진단 필드는 ANALYSIS_COHORT. 생검 필드는 BIOPSY_SITE. 발현은 쓰지 않았다. unmapped 진단은 at-risk에 넣지 않았다. 고유 장기 집합이 빈 부위(other)의 n_at_risk는 0이다.

{md(pd.read_csv(ROOT / "results/pog570_risk_by_site.tsv", sep="\\t"))}

## 14. 잠금 준수 확인
- `data/raw/POG570/` 안 파일은 POG570_README.txt, POG570_TPM_expression.txt.gz, Table_S1_Demographics.xlsx, Table_S2_Treatment.xlsx 네 개. 정책 HTML은 그 디렉터리 밖에 있다.
- 잠금 해시는 `config/locked_hashes.tsv`가 없을 때 한 번 썼고, 이후 검증했다.
- 발현 행렬을 진단 또는 생검 부위와 조인한 계산, 그림, 표는 없다.
- 발현 값을 읽은 곳은 유전자 ID 목록과 min/median/max 구조 요약뿐이다. 그 값은 진단과 붙지 않았다.
- unlock 없는 로더 호출은 예외로 끝났다.

## 15. 판단이 필요한 사항
1. MET500 `test==TRUE` 33개(평가 단위)를 포함할지 제외할지. 이번 실행은 포함.
2. MET500 SECR (평가 단위 20)을 unmapped로 둘지, 어느 장기로 둘지. 이번 실행은 unmapped.
3. MET500 MISC (평가 단위 39)를 unmapped로 둘지. 이번 실행은 unmapped.
4. MET500에 환자 ID가 없다. 부트스트랩 단위를 sample_source로 둘지. 이번 실행은 sample_source.
5. soft_tissue의 GTEx 참조가 Adipose - Subcutaneous와 Muscle - Skeletal 두 개일 때 비가중 평균을 쓸지. 이번 실행은 비가중 평균.
6. 시뮬레이션에서 Adipose - Subcutaneous의 고유 장기를 Sarcoma로 둘지. 이번 실행은 Sarcoma.
7. code 01 Recurrent Tumor TCGA-23-1023-01을 학습에서 뺄지. 이번 실행은 제외.
8. MET500 HNSC는 tissue가 여러 값이다. cohort로 장기를 정할지 tissue로 정할지. 이번 실행은 cohort.
9. POG570 진단 필드를 ANALYSIS_COHORT로 할지 TUMOUR_TYPE으로 할지. 이번 실행은 ANALYSIS_COHORT.
10. POG570 PANC 42에 췌장 신경내분비 종양이 들어 있다. Pancreas로 둘지 뺄지. 이번 실행은 PANC 전체를 Pancreas.
11. POG570 CNS-PNS와 BCC를 unmapped로 둘지. 이번 실행은 unmapped.
12. 생검 부위 사전은 초안이다. 사전에 없는 POG570 부위 48개는 other. other로 둘지 개별 부위로 나눌지.
13. SKCM 부위 필드를 `tumor_tissue_site`로 할지 `icd_o_3_site`로 할지. 이번 실행은 tumor_tissue_site.
14. 층화에 쓴 라이브러리 구분은 열 `sample_source`가 아니라 Sample_id의 capt/poly 토큰이다. 이 구분을 유지할지.
15. B0 유전자 5000개를 CV 밖에 한 번만 고를지, 폴드 안에서 다시 고를지. 이번 실행은 train 전체에서 한 번. B1 스케일러는 폴드 안에서 적합.
16. M1에서 남은 확률 합이 0이면 예측을 NA로 두고 오답으로 셀지. 이번 실행은 NA=오답.
17. Toil 원본 중앙값을 구간 0.001, 범위 [-15, 30] 히스토그램 중앙값으로 둘지 정확한 중앙값으로 둘지. 이번 실행은 히스토그램. 범위 밖 값은 0개. MET500 FPKM과 POG570은 정확한 중앙값.
18. README의 이용 조건(연구 목적, 재배포 금지, controlled dataset)과 웹 정책 페이지(CC BY 2.5, 상업 이용 포함) 중 어느 문장을 데이터 이용 조건으로 볼지.

## 16. 지시서와 다르게 한 것
- PyYAML 6.0.1을 가상환경에 설치해 YAML을 읽었다. numba 0.59.1로 빠른 KS를 계산했다. GEOparse는 쓰지 않았다.
- CV는 joblib 8프로세스 × BLAS 1스레드이다. 스크립트 시작 환경변수의 스레드 수는 16이고, 최종 적합은 16스레드이다. 8×16이 되지 않게 하려고 CV만 1스레드로 제한했다.
- phenotype은 latin-1로 읽었다.
- category, probemap, MET500 메타의 `.gz` URL이 403이라 비압축 파일을 저장했다.
- 정책 HTML은 `data/raw/POG570/`에 넣지 않고 `data/raw/pog570_data_release_policy.html`에 두어 잠금 디렉터리에는 공식 파일 4개만 있게 했다.
- 지시서의 묶음 행을 표준 부위 코드로 나눴다: colon_rectum, head_neck, bone과 bone_marrow, skin과 subcutaneous, peritoneum과 omentum. 고유 암종 집합은 묶음 행과 같다.
- GEO 발현은 받지 않았다. 헤더만 받았다.
- Toil 원본 중앙값은 히스토그램이다. 범위 밖 값은 0.
- GTEx `K-562` 세포주 70개는 기증자 ID를 만들지 않고 분할에서 제외했다.
- 빠른 KS는 검증을 통과해서 썼다. 원 함수로 전체 점수를 다시 계산하지는 않았다.

## 17. 산출물 목록
- `config/analysis_plan_A1.yaml`: §10 전에 잠근 분석 계획. 커밋 704890ea16d523911845ba17d9e43c14e70a9576.
- `config/split_tcga.tsv`, `config/split_gtex.tsv`, `config/met500_eval_samples.tsv`: 분할과 MET500 평가 단위.
- `config/mappings/*.tsv`: 대응표 초안.
- `config/locked_hashes.tsv`: POG570 원본 해시.
- `scripts/`: 다운로드, 감사, 분할, 행렬, KS, 측정, 로더.
- `logs/environment.txt`, `logs/disk_start.txt`, `logs/disk_after.txt`, `logs/download_log.tsv`, `logs/resource_usage.tsv`.
- `results/`: §10 표, KS 검증, 감사 빈도, 숙주 끌림 오류 목록.
- `data/processed/`: G, TPM 행렬, KS parquet. git에는 넣지 않는다.
- 보고서 본문 커밋: fb44fc642cb003549e9ae82f99dd52fc9e76f7c8. 이 해시를 본문에 적은 커밋이 그 다음 HEAD이다.
""")
    text = "\n".join(parts)
    dest = ROOT / "reports/A_보고서_1.md"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text)
    print("wrote", dest, "chars", len(text))


if __name__ == "__main__":
    main()
