# Complete annotated predictions — four historical datasets

These are byte-for-byte copies of the current report CSVs from `../outputs/<dataset>/`.
No rows or columns were removed, recalculated or renamed. No analysis was rerun.
These use the historical SRC/ERK2 cohorts, not the discarded re-curated datasets.

## Files and populations

| Dataset | Complete rows | Scored rows | Experimental positives (all / scored) | Columns |
|---|---:|---:|---:|---:|
| FGFR | 7,440 | 7,440 | 379 / 379 | 104 |
| HER2 | 5,870 | 5,870 | 114 / 114 | 107 |
| ERK2 | 537 | 509 | 16 / 8 | 115 |
| SRC | 873 | 873 | 101 / 101 | 115 |

Unlike `../ScoredVariants/`, this folder retains all 537 ERK2 rows, including
the 28 unscorable out-of-KD variants (8 experimental positives). Their raw model
scores and scenario calls are empty; they are not forced-negative predictions.

## Reading conventions and safeguards

- Experimental and predicted binary values use 1 for activating and 0 for non-activating.
- Numeric flags use 1/0; `sift_pp_keep` is stored as True/False. Prediction columns may show 0.0/1.0 because they permit missing values.
- Empty numeric prediction/score cells mean unavailable or scenario-ineligible, never automatically zero. An empty `exclusion_reason` instead means no annotation exclusion reason.
- Full-protein and KD residue coordinates are 1-based and bounds inclusive. Sequence strings use one-letter amino acids.
- `group_code` is an experimental label, not a family-feature code. For HER2 it is permissive; `evaluation_label` is conservative.
- GT HS is retrospective ground-truth-label assistance. It changes binary calls only, has last priority over zero overrides, and is not independent prospective model validation.
- SIFT/PolyPhen predict general sequence effects, not activation. Exclusion changes the evaluation population; AND/OR zero overrides keep all scorable rows and change calls only.
- No filter or hotspot changes the experimental labels or raw probabilities. Blank scenario predictions must not be interpreted as negatives.

## Models

- `GBDT`: investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2).
- `NGB`: NGBoost model.
- `SVM`: support-vector-machine classifier.
- `MNN`: five-checkpoint neural-network ensemble.
- `Classical_ensemble`: arithmetic mean of GBDT, NGB and SVM probabilities.
- `GBDT_MNN_ensemble`: arithmetic mean of investigator GBDT and MNN probabilities.

MNN uses frozen released weights and 50 BOS + 50 mutant-residue ESM PCA features.
Classical models additionally use the 67-column basic-feature tail (167 inputs total).
The PCA/feature matrices themselves are not contained in these CSVs.

## Scenario rules

| Scenario suffix | Definition |
|---|---|
| `Original` | All scorable rows; binary call is raw model probability >= 0.5. |
| `SIFT_PP` | Only scorable rows with sift_pp_keep=True; raw binary call. |
| `GT_HS` | All scorable rows; force the call to 1 at GT_HS_flag=1. |
| `SIFT_PP_GT_HS` | Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. |
| `Zero_AND` | All scorable rows; force the call to 0 at zero_AND_flag=1. |
| `Zero_AND_GT_HS` | All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. |
| `Zero_OR` | All scorable rows; force the call to 0 at zero_OR_flag=1. |
| `Zero_OR_GT_HS` | All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. |

All rules start with a raw-probability threshold of >=0.5. A hotspot never
restores an excluded or out-of-KD row to an evaluation population.

## Column-by-column dictionary

Every distinct column in the four copied files is listed explicitly below.
The availability column identifies files containing that field.

| Column | Available in | Meaning |
|---|---|---|
| `mutant_ID` | FGFR, HER2, ERK2, SRC | Normalized analysis identifier: protein symbol, underscore, amino-acid substitution; ERK2 uses MAPK1. |
| `protein` | FGFR, HER2, ERK2, SRC | Analysis protein/gene symbol (FGFR1–4, HER2, MAPK1 or SRC). |
| `mutation` | FGFR, HER2, ERK2, SRC | Substitution in one-letter notation, e.g. F19W: wild-type residue, 1-based full-protein position, mutant residue. |
| `group_code` | FGFR, HER2, ERK2, SRC | Original binary experimental input label, NOT a kinase-family encoding. Matches evaluation_label except HER2, where it is the permissive activation label. |
| `mutant_full_seq` | FGFR, HER2, ERK2, SRC | Mutated full-length protein sequence, retained as sequence/provenance metadata; the ESM input is the KD sequence, not this full protein. |
| `mutant_kd_seq` | FGFR, HER2, ERK2, SRC | Mutated selected kinase-domain sequence supplied to ESM for scorable rows. Out-of-KD ERK2 mutations leave this slice unchanged and are not scored. |
| `MNN_member_probability_sd` | FGFR, HER2, ERK2, SRC | Population SD (ddof=0) across the five sigmoid-transformed checkpoint logits. A descriptive member-spread measure, not a calibrated confidence interval. |
| `evaluation_label` | FGFR, HER2, ERK2, SRC | Unchanged primary experimental binary target: 1=activating, 0=non-activating under the supplied classification. HER2 uses conservative labels. Not a model prediction or hotspot override. |
| `inference_status` | FGFR, HER2, ERK2, SRC | scored or unscorable_outside_selected_KD. The latter does not imply an inactive prediction. |
| `position` | FGFR, HER2, ERK2, SRC | 1-based mutation position in the full protein, parsed from mutation. |
| `scorable` | FGFR, HER2, ERK2, SRC | 1 if a model score was generated for this variant; 0 for the 28 ERK2 mutations outside the agreed KD. |
| `GT_HS_flag` | FGFR, HER2, ERK2, SRC | Binary ground-truth hotspot flag at this protein/residue. A hotspot has >=50% activating tested substitutions at that residue. FGFR/HER2 use their pooled map (HER2 permissive labels), with training-map entries taking priority for both 0 and 1. ERK2/SRC use their own experimental-label maps, with the legacy unified map taking priority where keys overlap. Retrospective label-assisted metadata, not an independent prediction. |
| `GT_HS_origin` | FGFR, HER2, ERK2, SRC | Provenance label for the hotspot rule: retrospective_FGFR_HER2_plus_priority_training or retrospective_same_dataset_labels. The latter also follows legacy-map priority if a matching key exists. |
| `frozen_prior_HS_flag` | FGFR, HER2, ERK2, SRC | Hotspot flag from the legacy unified FGFR/HER2-plus-priority-training map alone, before supplementing with ERK2/SRC same-dataset maps. Not necessarily an independent external prior. |
| `in_original_training_variant_set` | FGFR, HER2, ERK2, SRC | 1 for an exact normalized mutant_ID match to an original training Kinase_mutation key; otherwise 0. This tests exact variants, not shared residue positions or all possible aliases. |
| `transcript_id` | FGFR, HER2, ERK2, SRC | Selected Ensembl transcript identifier used for sequence-matched variant annotation. |
| `protein_id` | FGFR, HER2, ERK2, SRC | Ensembl translated-protein identifier for the selected transcript; not the gene symbol. |
| `protein_hgvs` | FGFR, HER2, ERK2, SRC | Protein HGVS notation submitted for annotation, including protein version and three-letter amino-acid substitution. |
| `sift_prediction` | FGFR, HER2, ERK2, SRC | SIFT categorical effect call from Ensembl VEP, e.g. deleterious or tolerated (including low-confidence categories if present). This does not directly predict kinase activation. |
| `sift_score` | FGFR, HER2, ERK2, SRC | Raw SIFT score on 0–1 scale; lower means more likely deleterious. Stored as returned, not an activation probability; filtering uses categorical calls. |
| `polyphen_prediction` | FGFR, HER2, ERK2, SRC | PolyPhen-2 categorical effect call from Ensembl VEP, e.g. benign, possibly_damaging or probably_damaging. Not an activation label. |
| `polyphen_score` | FGFR, HER2, ERK2, SRC | Raw PolyPhen-2 score on 0–1 scale; higher means more likely damaging. Not an activation probability; filtering uses categorical calls. |
| `annotation_status` | FGFR, HER2, ERK2, SRC | OK when sequence/substitution-matched annotation is accepted; ERR-prefixed values denote unavailable, mismatched or conflicting annotation. All four current files have OK annotations. |
| `sift_pp_keep` | FGFR, HER2, ERK2, SRC | Boolean annotation-only exclusion-rule retention: annotation_status=OK AND SIFT is not tolerated/tolerated_low_confidence AND PolyPhen is not benign. Scenario eligibility additionally requires scorable=1. |
| `exclusion_reason` | FGFR, HER2, ERK2, SRC | Reason for failing the annotation-only exclusion rule; SIFT_tolerated takes reporting precedence over PolyPhen_benign if both apply. Empty means no exclusion reason; ERR reasons denote annotation problems if present. Does not indicate exclusion from Original scoring. |
| `annotation_origin` | FGFR, HER2, ERK2, SRC | Annotation provenance; these four files use fresh_Ensembl_VEP, not paper-supplied annotations. |
| `retrieved_utc` | FGFR, HER2, ERK2, SRC | UTC retrieval timestamp recorded by the annotation run, not the date of this copy. |
| `Original_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to Original evaluation, otherwise 0. All scorable rows; binary call is raw model probability >= 0.5. |
| `SIFT_PP_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to SIFT_PP evaluation, otherwise 0. Only scorable rows with sift_pp_keep=True; raw binary call. |
| `GT_HS_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to GT_HS evaluation, otherwise 0. All scorable rows; force the call to 1 at GT_HS_flag=1. |
| `SIFT_PP_GT_HS_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to SIFT_PP_GT_HS evaluation, otherwise 0. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. |
| `GBDT_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2), on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `GBDT_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2): 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `NGB_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from NGBoost model, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `NGB_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for NGBoost model: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `SVM_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from support-vector-machine classifier, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `SVM_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for support-vector-machine classifier: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `MNN_member_1_logit` | FGFR, HER2, ERK2, SRC | Raw unbounded pre-sigmoid output of neural checkpoint 1 (of five), not a probability. |
| `MNN_member_2_logit` | FGFR, HER2, ERK2, SRC | Raw unbounded pre-sigmoid output of neural checkpoint 2 (of five), not a probability. |
| `MNN_member_3_logit` | FGFR, HER2, ERK2, SRC | Raw unbounded pre-sigmoid output of neural checkpoint 3 (of five), not a probability. |
| `MNN_member_4_logit` | FGFR, HER2, ERK2, SRC | Raw unbounded pre-sigmoid output of neural checkpoint 4 (of five), not a probability. |
| `MNN_member_5_logit` | FGFR, HER2, ERK2, SRC | Raw unbounded pre-sigmoid output of neural checkpoint 5 (of five), not a probability. |
| `MNN_logit` | FGFR, HER2, ERK2, SRC | Arithmetic mean of the five checkpoint logits, before sigmoid. |
| `MNN_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from five-checkpoint neural-network ensemble, on 0–1 scale. Computed as sigmoid(mean of five member logits), not mean of member probabilities. Never replaced by hotspot/zero overrides. |
| `MNN_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for five-checkpoint neural-network ensemble: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `Classical_ensemble_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from arithmetic mean of GBDT, NGB and SVM probabilities, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `GBDT_MNN_ensemble_probability` | FGFR, HER2, ERK2, SRC | Unmodified activating-class probability from arithmetic mean of investigator GBDT and MNN probabilities, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `Classical_ensemble_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for arithmetic mean of GBDT, NGB and SVM probabilities: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `GBDT_MNN_ensemble_prediction` | FGFR, HER2, ERK2, SRC | Baseline binary call for arithmetic mean of investigator GBDT and MNN probabilities: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `zero_AND_flag` | FGFR, HER2, ERK2, SRC | 1 if annotation_status=OK AND SIFT is tolerated/tolerated_low_confidence AND PolyPhen is benign; flags a forced-negative call in Zero_AND scenarios, not row removal. |
| `zero_OR_flag` | FGFR, HER2, ERK2, SRC | 1 if annotation_status=OK AND either SIFT is tolerated/tolerated_low_confidence OR PolyPhen is benign; flags a forced-negative call in Zero_OR scenarios, not row removal. |
| `annotation_missing_or_conflicting` | FGFR, HER2, ERK2, SRC | 1 when annotation_status is not OK; otherwise 0. Such rows are retained by Original and full-population zero-override rules but not the exclusion scenario (if scorable). |
| `GBDT_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Original_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_SIFT_PP_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_AND_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to Zero_AND evaluation, otherwise 0. All scorable rows; force the call to 0 at zero_AND_flag=1. |
| `GBDT_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_AND_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_AND_GT_HS_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to Zero_AND_GT_HS evaluation, otherwise 0. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. |
| `GBDT_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_OR_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to Zero_OR evaluation, otherwise 0. All scorable rows; force the call to 0 at zero_OR_flag=1. |
| `GBDT_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_OR_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_OR_GT_HS_eligible` | FGFR, HER2, ERK2, SRC | 1 if row contributes to Zero_OR_GT_HS evaluation, otherwise 0. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. |
| `GBDT_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for investigator-trained gradient-boosted decision-tree model (not collaborator GBDT-2) under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for NGBoost model under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for support-vector-machine classifier under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for five-checkpoint neural-network ensemble under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2, SRC | Binary call for arithmetic mean of investigator GBDT and MNN probabilities under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `activation_conservative` | HER2 | HER2 binary experimental GOF label requiring concordant activity calls in both technical replicates; primary endpoint. |
| `activation_permissive` | HER2 | HER2 broader supplied binary experimental GOF label, including single-replicate calls; equals group_code. |
| `permissive_only_ambiguous` | HER2 | 1 when HER2 is permissive-positive but conservative-negative; otherwise 0. These rows remain negatives in the primary evaluation. |
| `source_ID` | ERK2, SRC | Original supplied spreadsheet identifier, e.g. ERK2_F19W, before normalization to MAPK1. |
| `source_row_1_based` | ERK2, SRC | Original Excel row number including the header (first data row is 2). |
| `fidelity` | ERK2, SRC | Fidelity metadata copied as an integer from the original spreadsheet Fidelity column; not a model confidence measure and not re-derived here. |
| `kd_start` | ERK2, SRC | Inclusive 1-based start of the selected inference KD: ERK2 19, SRC 270. |
| `kd_end` | ERK2, SRC | Inclusive 1-based end of the selected inference KD: ERK2 353, SRC 523. |
| `kd_position_1_based` | ERK2, SRC | Full-protein position minus kd_start plus 1. For unscorable ERK2 rows this can be <=0 or beyond KD length; it is not a valid scored index in those rows. |
| `core_kd_start` | ERK2, SRC | Inclusive 1-based internal database core-KD start: ERK2 25, SRC 270. |
| `core_kd_end` | ERK2, SRC | Inclusive 1-based internal database core-KD end: ERK2 313, SRC 523. |
| `in_internal_core` | ERK2, SRC | 1 if position lies within core_kd_start..core_kd_end, otherwise 0; distinct from membership of the selected extended ERK2 inference KD. |
| `domain_definition` | ERK2, SRC | Selected-domain policy: supplied_extended_kinase_domain for ERK2; internal_annotated_kinase_domain for SRC. |
| `uniprot_accession` | ERK2, SRC | Sequence reference accession: P28482 for MAPK1/ERK2, P12931 for SRC. |

## Provenance and integrity

Documentation was checked against `../evaluate.py`, `../../evaluate_analysis.py`,
`../../run_inference.py`, `../../prepare_analysis.py` and `../../annotate_vep.py`.
The original spreadsheet definition of `Fidelity` is not reinterpreted as confidence.
Files were copied without reserialization and their SHA-256 hashes checked:

- `FGFR_complete_annotated_predictions.csv`: `24ad7266377d904aa664a8b5cf7cd6ae58a76b5f1f0e99362727b5f42e76dd90`
- `HER2_complete_annotated_predictions.csv`: `f95febba653459c2dbe847e844745c811e85c4bdea71de19e52bfdc150e2b1c1`
- `ERK2_complete_annotated_predictions.csv`: `dbd70f08888ca6c6d6a09277bfc0a0f8304eaa67c9b13669afcdaa46e3cf4469`
- `SRC_complete_annotated_predictions.csv`: `09135bf5c656476c640a0eead1d22aff12de473eeafb36ea063c971ee5b21836`

Data dictionary coverage: 118 of 118 distinct columns documented.
