# A7 잠금 해제 덧붙임

잠금 파일 `config/unlock_A7.md`는 고치지 않는다. 이 파일은 A 지시서 8 §1의 6번 결정에 따라 시간순을 적는다.

| 커밋 | 시각 | 내용 |
|---|---|---|
| `84086a1a` | 2026-10-02 12:26:30 +0900 | 코드 동결 |
| `198007aa` | 2026-10-02 12:27:10 +0900 | 출력 경로 (잠금 파일에 적힌 코드 커밋) |
| `edd32463` | 2026-10-02 12:27:53 +0900 | 잠금 해제 |
| — | — | 첫 `--real`. 빈 `standard_site`(NaN)로 `score_base`에서 종료. 가설 표 없음 |
| `a7272492` | 2026-10-02 12:29:49 +0900 | 빈 값 처리 수정, `frozen_code_A7.tsv` 해시 갱신 |
| — | — | 두 번째 `--real` (사전등록 결과) |
| `8f788feb` | 2026-10-02 12:34:15 +0900 | 출력 해시 |

실행에 쓴 `scripts/24_aux_confirm.py` SHA-256 `9d1fc69218eca2c37048c66a6600d039f155f9da783238dcb1c9ae21eedb2806`.

`git diff 198007aa a7272492 --stat`:

```
 config/frozen_code_A7.tsv |  2 +-
 config/unlock_A7.md       | 16 ++++++++++++++++
 scripts/24_aux_confirm.py |  5 ++++-
 3 files changed, 21 insertions(+), 2 deletions(-)
```

이 범위에는 잠금 해제 커밋 `edd32463`이 들어 있어 `config/unlock_A7.md` 16줄 추가가 함께 나온다. 그 파일을 빼면 바뀐 것은 `score_base`의 빈 값 처리와 `frozen_code_A7.tsv`의 `24_aux_confirm.py` 한 줄이다.

`scripts/24_aux_confirm.py` 바뀐 줄:

```
-    labels["library"] = labels["library"].fillna("")
+    # Empty label fields are NaN under read_csv even with dtype=str. Excluded rows have an empty site.
+    for column in ("library", "standard_site", "organ", "native_organs", "patient_id", "sample_id"):
+        if column in labels.columns:
+            labels[column] = labels[column].fillna("").astype(str)
```

`config/frozen_code_A7.tsv` 바뀐 줄:

```
-scripts/24_aux_confirm.py	code	54839	b4b63c95a042531f0f195304c731c6fd18635557fb2e798967faa9dbe39b5664
+scripts/24_aux_confirm.py	code	55094	9d1fc69218eca2c37048c66a6600d039f155f9da783238dcb1c9ae21eedb2806
```

이 처리는 A 지시서 7 §3.4(표가 나오기 전 실패는 고쳐 다시 실행)에 해당한다.
