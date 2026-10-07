# A 사전등록 3 — POG570 확인 검정

이 문서는 POG570 발현으로 특징을 계산하기 전에 고정한다. 정답 장기와 예측의 결합은 이 문서의 커밋 뒤 다음 지시서에서 한 번만 한다.

## 이전 단계 결정

보고서 2 §10의 번호.

| 번호 | 결정 |
| --- | --- |
| 1 | HPA 중복 Gene name은 nTPM을 합한 뒤 log2(sum+1) |
| 2 | IF 제거 개수는 floor(n×q) |
| 3 | SA 숙주 추출은 조직 균등 후 그 조직의 샘플 균등 |
| 4 | V0는 1단계 규칙. proxy 조직을 넣지 않음 |
| 5 | NC에서 종양 클래스 확률 합이 0이면 예측 NA, 오답 |
| 6 | conformal에서 k > n_cal이면 임계값 +inf |
| 7 | 새 방법 동률은 MET500 top-1, 그다음 host_rate, 그다음 method_id 알파벳순 |
| 8 | β 3분위는 안정 정렬 140/140/139. POG570 구간은 MET500 경계로 고정: β = 0, 0 < β < 0.0327, β ≥ 0.0327. 0.0327은 MET500 T3 최솟값 0.032692468741067175를 소수 넷째 자리에서 반올림한 값 |
| 9 | soft_tissue 평균은 LD·증강이 샘플 풀 평균 후 합 1e6. V0는 두 조직 평균 벡터의 비가중 평균 |

## 후보

- 주 비교: SA-Z 대 BASE-Z. 2단계 규칙으로 MET500 top-1이 가장 높은 새 방법이다. top-1은 0.7276887871853547이었다.
- 보조 비교. 다중성 보정 없음.
  - SC-Z, LD-Z, NC-Z, M1-Z 대 BASE-Z
  - SA-K, V0-K 대 BASE-K
  - BASE-K 대 BASE-Z
- 모델·참조·대응표 SHA-256은 `config/frozen_models_A3.tsv`.

| path | role | size_bytes | sha256 |
| --- | --- | --- | --- |
| results/stage2/models/BASE_Z.joblib | model | 1281982 | 9392f6252854b91f37aad20a3114200d951d8585c284f809d8d77ddaad776b32 |
| results/stage2/models/BASE_K.joblib | model | 2317125 | 4a94c8ec883870ccd6991d82e795144cb598e83b1d62438b083b24a64406d2d8 |
| results/stage2/models/SA_Z.joblib | model | 1402359 | 7c39219d35349d83792b012f914f7c2ab237b83f3c3956accafe8781bde6ddea |
| results/stage2/models/SA_K.joblib | model | 2317103 | c5575dbee2b93e799351beeb0074d62f58c5f239dd713f62945f4c26547e8ee4 |
| results/stage2/models/NC_Z.joblib | model | 1802884 | 5aceaad78b52437e810a3fd1ae61da0e29e3bf0ff7ccc5e3e8dc5d243120880b |
| results/stage2/models/SC_adrenal_Z.joblib | model | 1402376 | a2415c7defbab04d7dc9cf38a60f277d508d5d69cc25ae2f1c2504a10689d8a6 |
| results/stage2/models/SC_bone_marrow_Z.joblib | model | 1402380 | 2fabbdd61ab1d0d1668ce09e21aa9a38471a98323432a825ffe74b6c4d2c95af |
| results/stage2/models/SC_brain_Z.joblib | model | 1402374 | 3aa0930e3ec40fc8f011ad2e0a408f83ff5aaf028e682d08924060fd5e6114f4 |
| results/stage2/models/SC_liver_Z.joblib | model | 1402374 | e64214634dfeb622a927ac446d7aa99e251ac52017ff36d716464be6f70ff804 |
| results/stage2/models/SC_lung_Z.joblib | model | 1402373 | 825a88b37986806cf1989130a87c05056091f4fa607fb817a4c2bd5a1a299bad |
| results/stage2/models/SC_lymph_node_Z.joblib | model | 1402379 | da0c29257c067dc48437cb864f95f0abfcf63a643f7a1b2dd7239e106777bc95 |
| results/stage2/models/SC_omentum_Z.joblib | model | 1402376 | 99efa9d4e3c28456c35f08964616c95c85ae1e02f31cba1b76fc8cfa69524f4b |
| results/stage2/models/SC_skin_Z.joblib | model | 1402373 | d344bcfa5d9ec7bd3a710be2f48af00b0473858838f2a553540389ec267fe62d |
| results/stage2/models/SC_soft_tissue_Z.joblib | model | 1402380 | 32d47d80c9bf53271b8e71b29d02222daba8fd058b2135777bde49670ae4b850 |
| config/mappings/met500_cohort_to_organ.tsv | mapping | 292 | 8c6260d8415c77c1a420bdd5ef8e03ed72b58fed637fdc2b6ccdb7746a3d3804 |
| config/mappings/pog570_cohort_to_organ.tsv | mapping | 365 | a3da27c3e5ca67e5675e07c6c920fe80367e242fcf0c5895ff5fe3c47b5dfded |
| config/mappings/project_to_organ.tsv | mapping | 425 | 48db237b3b50416891bdb8ffb3fdfa6ecf76f0ef1d01b195bf92a80a47523cf7 |
| config/mappings/site_dictionary.tsv | mapping | 2702 | b688f437a4abcbdf2cf6faac22613a7f72f9fb4517a4f27866603bdd1c817f77 |
| config/mappings/site_native.tsv | mapping | 1081 | ad38d55db462e1cde354c2e941a25fbf151fbd89bd5c927c043a601b205387ab |
| config/mappings/tcga_disease_to_project.tsv | mapping | 1162 | 8323787f9cef7b248aa7f5fe6f7af6af3dccc6516e80936dcc79a105711c98fa |
| results/B0_genes.txt | b0_genes | 31057 | 51aea84ef00fac6d7daeb6e0713d870367071f85280245bd847e28706703b406 |
| data/processed/genes/G_symbols.txt | G | 111844 | b9085bf65f64c6140a2b5e32cf4e2074f06b810f87ccef4777224538fd5aac09 |
| data/processed/genes/gene_sets.npz | ks_sets | 3189791 | dd0104c81f0794545730d93d68f487c434d26e6518f42e44bfaff22ac871168e |
| results/stage3/references.npz | ld_v0_references | 3712241 | 30a343d76082cc18669577d24308880edd602f31b4a87071606a75d133cdd5a1 |

## 평가 집합

라벨 파일 `config/pog570_eval_labels.tsv`. SHA-256 `5897cb0db827ca6f3d95504a6e44c84ed6f929ff5a2430c97d882a308dd6f4c2`.

예외 규칙을 코호트 표보다 먼저, 1번부터 7번까지 적용한다. 처음 맞는 규칙만 남긴다. 맞는 예외가 없으면 코호트 표(rule = basic)를 쓴다.

1. TUMOUR_TYPE 또는 HISTOLOGICAL_TYPE에 `neuroendocrine` 또는 `carcinoid`가 있으면 제외한다. 대소문자는 무시한다. `small cell`도 같은 방식이다. 다만 그 일치의 바로 앞이 `non-` 또는 `non `이면 제외하지 않는다. `non-small cell`은 이 규칙에 걸리지 않는다. `net`은 검색하지 않는다. HISTOLOGICAL_TYPE을 strip·대문자로 만든 값이 LUNE, CSCLC, SCLC, ALUCA, LNET, PANET, SCCO이면 제외한다. 코드는 부분 문자열이 아니다.
2. ANALYSIS_COHORT가 LUNG이고 HISTOLOGICAL_TYPE이 PEMESO, PLMESO, PLEMESO이면 Mesothelioma.
3. CNS-PNS: GBM, AODG, ODG, PAST, GNG → Brain. MPNST → Sarcoma. PGNG → Adrenal. MNG, GRCT, APE, MPE, EPM, SFTCNS, PNET → 제외. 이 목록에 없는 코드는 NA, rule = 3-unlisted.
4. MISC: ESCC → Esophagus. NSGCT → Testis. 그 외 코드는 제외.
5. HNSC: OPHSC → HeadNeck. HNMUCM → Melanoma. NPC, ODGC, CSCC, SNA → 제외. 그 외 HNSC는 코호트 표로 넘어간다.
6. SARC: ATRT, MYEC → 제외. 그 외는 Sarcoma.
7. OV: OCS → 제외. 그 외는 Ovary. SCCO는 1번에서 이미 제외된다.
- 코호트 표: BRCA Breast, COLO Colorectal, LUNG Lung, SARC Sarcoma, PANC Pancreas, OV Ovary, CHOL Biliary, SKCM Melanoma, STAD Stomach, UCEC Uterus, LYMP Lymphoid, ESCA Esophagus, HNSC HeadNeck, UVM UvealMelanoma, KDNY Kidney, ACC Adrenal, CERV Cervix, THCA Thyroid, THYM Thymus, PRAD Prostate, HCC Liver, BLCA Bladder. SECR와 BCC는 제외. CNS-PNS와 MISC는 이 표에 없다.

생검 부위는 1단계 `site_dictionary.tsv`의 pog570_biopsy_site. 고유 장기 집합은 1단계 `site_native.tsv`.

집합 정의.

- 평가 집합: organ이 exclude도 NA도 아닌 샘플. n = 512.
- 위험 집합: 평가 집합 중 표준 부위의 고유 장기 집합이 비어 있지 않고 정답 장기가 그 집합 밖인 샘플. n = 378.
- 고유 장기 정답 집합: 정답 장기가 그 집합 안인 샘플. n = 91.
- 제외 58. NA 0. 평가 집합 중 고유 장기 집합이 빈 샘플은 43이고, standard_site가 other인 평가 샘플 43과 같다.

TUMOUR_CONTENT 3분위는 평가 집합 512개의 값에 `qcut` 3, duplicates=drop을 적용한 구간이다. 닫힘은 right. 확인 검정은 이 경계를 다시 계산하지 않는다.

| tertile | left | right | closed | n |
| --- | --- | --- | --- | --- |
| T1 | 16.999 | 48.0 | right | 173 |
| T2 | 48.0 | 66.667 | right | 168 |
| T3 | 66.667 | 100.0 | right | 171 |

## 가설

고정 순서. 각 단측 α = 0.05. 환자당 샘플은 1개다. 부트스트랩은 그 가설의 분석 집합을 환자 단위로 2,000회 재추출한다. 시드 20261001. 백분위 2.5와 97.5. Generator 호출 순서는 H1, H2, H3이다. 차이는 SA-Z 지표 − BASE-Z 지표이다.

| 순서 | 가설 | 검정 |
| --- | --- | --- |
| H1 | 위험 집합에서 SA-Z의 숙주 끌림 오류율이 BASE-Z보다 낮다 | 숙주 끌림 오류 지시자의 정확 McNemar 단측. 유리한 불일치 수는 BASE만 숙주 오류인 수. p는 Binomial(n_disc, 0.5)에서 그 수가 관측값 이상일 확률. n_disc = 0이면 p = 1 |
| H2 | 평가 집합에서 SA-Z의 장기 top-1이 BASE-Z보다 높다 | 정확 McNemar 단측. 유리한 불일치 수는 SA만 정답인 수. H1이 유의할 때만 검정한다 |
| H3 | 고유 장기 정답 집합에서 SA-Z의 top-1이 BASE-Z보다 10%p 넘게 낮지 않다 | 같은 H3 부트스트랩 차이의 5백분위수가 −0.10보다 크다. H2가 유의할 때만 판정한다 |

H1·H2·H3 모두 차이와 양측 95% 신뢰구간을 보고한다. H3은 단측 5백분위수도 보고한다.

예측 NA는 오답이다.

## 보조 분석

다중성 보정 없음. 점추정만 보고한다.

- 방법별 평가 집합 지표: top-1, top-3, macro-F1, host_rate, 고유 장기 정답 top-1. 비교 목록은 위의 보조 비교와 두 BASE이다.
- 부위 묶음: liver, lymph_node, lung, soft_tissue, remainder. remainder는 그 넷이 아닌 표준 부위이다. 각 방법의 n, top-1, n_at_risk, host_rate.
- METASTATIC_OR_RECURRENCE가 Metastatic인 평가 샘플만. 같은 전체 지표.
- TUMOUR_CONTENT 구간 T1, T2, T3. 위의 경계. 각 방법의 n, top-1, n_at_risk, host_rate.
- LD β 구간: 0, (0, 0.0327), ≥ 0.0327. 유한하지 않으면 NA. 음수이면 negative. BASE-Z와 SA-Z의 n, top-1, n_at_risk, host_rate.
- conformal. 임계값은 MET500 평가 437개 전체로 정한다. 방법 BASE-Z, SA-Z. 변형 global과 Mondrian {liver, lymph_node, other}. α = 0.1, 0.2. 점수 s = 1 − p̂(정답 장기). k = ceil((n_cal+1)×(1−α)). k > n_cal이면 임계값 +inf. 집합은 1 − p(o) ≤ 임계값인 장기. 지표는 포함률, 평균 집합 크기, 단일 원소 비율, 단일 원소 정확도. 빈 집합의 크기는 0이다.

사전등록에 없는 분석은 사후 분석으로 따로 보고한다.

## 처리

- POG570 TPM은 버전 없는 Ensembl을 1단계와 같은 HGNC 기호로 바꾼다. 한 기호에 여러 Ensembl이 오면 TPM을 합한다. G에 없는 기호는 버리고, G에 없는 행은 0으로 채운 뒤 합 1e6으로 맞춘다.
- Z는 G 전체 순위-정규 뒤 1단계 B0 유전자 5,000개. K는 1단계 KS 구현.
- 부위 모델이 없는 부위의 SC는 SA를 쓴다. LD, M1, V0는 1·2단계 규칙을 따른다.
- 예측 파일은 570명 전체이다. 평가 집합 필터는 라벨 파일을 읽을 때만 적용한다.
- 스크립트의 OMP, OPENBLAS, MKL, NUMEXPR, NUMBA 스레드는 16이다.
- 가짜 라벨 시험의 정답 장기는 시드 1이다. 분포는 보고서 1의 ANALYSIS_COHORT 단독 빈도를 코호트 표의 장기로 바꾼 다항분포이다. CNS-PNS, MISC, SECR, BCC는 그 분포에 넣지 않는다. 가중 합은 525이다. 추출 순서는 Table S1 행 순서이다. 샘플별 정답 장기는 그 추출의 입력이 아니다.

## 검정력

정확 McNemar 단측 α = 0.05. 시뮬레이션 10,000회. 시드 20261001. 시나리오 순서는 H1 전체, H1 절반, H2 전체, H2 절반이다. 각 샘플은 확률 p_favor로 가설에 유리한 불일치, p_against로 불리한 불일치, 그 밖으로 일치를 낸다. 절반 효과는 p_favor만 절반이고 p_against는 그대로이다.

| 시나리오 | 가설 | n | p_favor | p_against | power | n_significant |
| --- | --- | --- | --- | --- | --- | --- |
| MET500 비율 | H1 | 378 | 34/361 = 0.09418282548476455 | 1/361 = 0.002770083102493075 | 1.0 | 10000 |
| 효과 절반 | H1 | 378 | 17/361 = 0.04709141274238227 | 1/361 = 0.002770083102493075 | 0.9959 | 9959 |
| MET500 비율 | H2 | 512 | 29/437 = 0.06636155606407322 | 12/437 = 0.02745995423340961 | 0.8727 | 8727 |
| 효과 절반 | H2 | 512 | 29/874 = 0.03318077803203661 | 12/437 = 0.02745995423340961 | 0.0985 | 985 |

## 이전 단계 변경 이력

- 1단계 계획 커밋 704890ea16d523911845ba17d9e43c14e70a9576 뒤에 분석 설정 코드는 바뀌지 않았다. 그 뒤 커밋은 보고서와 보고서 해시 기록이다.
- 2단계 계획 커밋 54b81c045a117cb19d8d4bb9721a81bc1e6b34f1 뒤에 `scripts/09_eval.py`의 β 3분위가 바뀌었다. 커밋된 호출은 `pd.qcut(..., 3, labels=[T1,T2,T3], duplicates="drop")`이었다. 유한 β 419개 중 167개가 0이라 ValueError가 났다. 커밋 5e8f0c3a2b54d086f3019642a91397637b87238b는 안정 정렬로 140, 140, 139를 나눈다. 계획 YAML은 고치지 않았다. 09_eval을 처음부터 다시 실행했다. `logs/time_09.txt`는 재실행 기록이다.
- HPA 중복 심볼의 nTPM 합산은 2단계 계획 커밋 전에 다시 계산했다.
