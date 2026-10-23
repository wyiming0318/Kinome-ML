# Independent-validation datasets and predictions

This directory contains **six CSV files for three datasets: FGFR, HER2 and ERK2/MAPK1**.
Each dataset has a compact scored-variant table and a complete annotated prediction table.
SRC data are not included in this directory.

## Files and populations

| Dataset | Scored variants: rows / columns | Complete predictions: rows / columns | Experimental positives: scored / complete |
|---|---|---|---|
| FGFR | [FGFR_scored_variants.csv](FGFR_scored_variants.csv): 7,440 / 8 | [FGFR_complete_annotated_predictions.csv](FGFR_complete_annotated_predictions.csv): 7,440 / 104 | 379 / 379 |
| HER2 | [HER2_scored_variants.csv](HER2_scored_variants.csv): 5,870 / 11 | [HER2_complete_annotated_predictions.csv](HER2_complete_annotated_predictions.csv): 5,870 / 107 | 114 / 114 |
| ERK2 | [ERK2_scored_variants.csv](ERK2_scored_variants.csv): 509 / 19 | [ERK2_complete_annotated_predictions.csv](ERK2_complete_annotated_predictions.csv): 537 / 115 | 8 / 16 |

### Which files should I use?

- **`*_scored_variants.csv`** lists exactly the variants with model scores, before SIFT/PolyPhen subset exclusion. It includes experimental labels, identifiers, mutation positions, mutant full-protein/KD sequences and available source/domain metadata. It does **not** contain model predictions, SIFT/PolyPhen annotations or hotspot/scenario outputs.
- **`*_complete_annotated_predictions.csv`** includes the full retained population, those metadata fields, raw model probabilities and baseline calls, five neural-checkpoint logits, SIFT/PolyPhen annotations, hotspot flags, and eligibility/calls for all eight evaluation scenarios.
- The scored files are compact projections of the scorable rows in the corresponding complete files, not additional cohorts. Do not concatenate each pair as independent observations.
- Neither file type contains the 1,280-dimensional ESM vectors, PCA input matrices, 67-column model basic-feature vectors, model weights or aggregate performance tables.

**ERK2:** the scored file contains 509 variants (8 positives). The complete file retains all
537 source variants (16 positives), including 28 outside the selected KD, residues 19–353.
Those 28 include 8 positives; their model scores and scenario predictions are empty, not
forced-negative calls. They are absent from the scored file. The complete files for FGFR
and HER2 have no additional unscorable rows.

## Experimental labels and conventions

- `evaluation_label` is the primary experimental target: 1 = activating, 0 = non-activating under the supplied binary classification. It is not a model prediction.
- HER2 uses `activation_conservative` as its primary endpoint (114 positives). The broader `activation_permissive` label has 438 positives; the 324 permissive-only variants remain negatives under the conservative endpoint. Both labels are supplied in both HER2 files.
- `group_code` matches the experimental target for FGFR/ERK2, but represents the permissive label for HER2. It is not a kinase-family encoding.
- ERK2 uses normalized `MAPK1` identifiers in `mutant_ID`; the original `ERK2` identifiers remain in `source_ID`.
- Residue positions are one-based and KD bounds are inclusive. Sequence strings use one-letter amino-acid codes.
- Flags use 1/0; `sift_pp_keep` is stored as True/False. Prediction columns can appear as 0.0/1.0 to allow missing values.
- Empty score/prediction cells mean unscorable or scenario-ineligible, not zero. Empty `exclusion_reason` instead means no annotation exclusion reason.

## Models and interpretation

- `GBDT`, `NGB` and `SVM`: frozen gradient-boosted decision-tree, NGBoost and support-vector-machine models.
- `MNN`: the five-checkpoint memory-based neural-network ensemble. Its probability is the sigmoid of the mean member logit, not the mean member probability.
- `Classical_ensemble`: arithmetic mean of GBDT, NGB and SVM probabilities.
- `GBDT_MNN_ensemble`: arithmetic mean of GBDT and MNN probabilities.

The released MNN uses 50 BOS + 50 mutation-residue ESM PCA features (100 inputs).
Classical models additionally use 67 family/basic features (167 inputs). These input
vectors and saved preprocessing are not distributed in this folder.

The `Original` scenario is the unmodified external model evaluation. **GT HS scenarios
are retrospective, label-assisted analyses**, because the hotspot maps use experimental
labels from the evaluation datasets and the historical training prior. They should not
be interpreted as independent prospective validation. SIFT/PolyPhen calls describe general
sequence effects, not kinase activation or experimental folding status. Exclusion changes
the evaluated population; AND/OR zero overrides change binary calls without removing
scorable variants. Hotspot overrides have last priority. No scenario changes experimental
labels or raw model probabilities.

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

The original column-explanation table is retained below and updated for the datasets
actually present. **“Available in” refers to the complete annotated prediction files.**
The scored files contain only the compact metadata subset described above; every column
in those files is also explained by this table. There are 118 distinct columns across the
six CSVs.

| Column | Available in | Meaning |
|---|---|---|
| `mutant_ID` | FGFR, HER2, ERK2 | Normalized analysis identifier: protein symbol, underscore, amino-acid substitution; ERK2 uses MAPK1. |
| `protein` | FGFR, HER2, ERK2 | Analysis protein/gene symbol (FGFR1–4, HER2 or MAPK1). |
| `mutation` | FGFR, HER2, ERK2 | Substitution in one-letter notation, e.g. F19W: wild-type residue, 1-based full-protein position, mutant residue. |
| `group_code` | FGFR, HER2, ERK2 | Original binary experimental input label, NOT a kinase-family encoding. Matches evaluation_label except HER2, where it is the permissive activation label. |
| `mutant_full_seq` | FGFR, HER2, ERK2 | Mutated full-length protein sequence, retained as sequence/provenance metadata; the ESM input is the KD sequence, not this full protein. |
| `mutant_kd_seq` | FGFR, HER2, ERK2 | Mutated selected kinase-domain sequence supplied to ESM for scorable rows. Out-of-KD ERK2 mutations leave this slice unchanged and are not scored. |
| `MNN_member_probability_sd` | FGFR, HER2, ERK2 | Population SD (ddof=0) across the five sigmoid-transformed checkpoint logits. A descriptive member-spread measure, not a calibrated confidence interval. |
| `evaluation_label` | FGFR, HER2, ERK2 | Unchanged primary experimental binary target: 1=activating, 0=non-activating under the supplied classification. HER2 uses conservative labels. Not a model prediction or hotspot override. |
| `inference_status` | FGFR, HER2, ERK2 | scored or unscorable_outside_selected_KD. The latter does not imply an inactive prediction. |
| `position` | FGFR, HER2, ERK2 | 1-based mutation position in the full protein, parsed from mutation. |
| `scorable` | FGFR, HER2, ERK2 | 1 if a model score was generated for this variant; 0 for the 28 ERK2 mutations outside the agreed KD. |
| `GT_HS_flag` | FGFR, HER2, ERK2 | Binary ground-truth hotspot flag at this protein/residue. A hotspot has >=50% activating tested substitutions at that residue. FGFR/HER2 use their pooled map (HER2 permissive labels), with training-map entries taking priority for both 0 and 1. ERK2 uses its own experimental-label map, with the legacy unified map taking priority where keys overlap. Retrospective label-assisted metadata, not an independent prediction. |
| `GT_HS_origin` | FGFR, HER2, ERK2 | Provenance label for the hotspot rule: retrospective_FGFR_HER2_plus_priority_training or retrospective_same_dataset_labels. The latter also follows legacy-map priority if a matching key exists. |
| `frozen_prior_HS_flag` | FGFR, HER2, ERK2 | Hotspot flag from the legacy unified FGFR/HER2-plus-priority-training map alone, before supplementing with the ERK2 same-dataset map. Not necessarily an independent external prior. |
| `in_original_training_variant_set` | FGFR, HER2, ERK2 | 1 for an exact normalized mutant_ID match to an original training Kinase_mutation key; otherwise 0. This tests exact variants, not shared residue positions or all possible aliases. |
| `transcript_id` | FGFR, HER2, ERK2 | Selected Ensembl transcript identifier used for sequence-matched variant annotation. |
| `protein_id` | FGFR, HER2, ERK2 | Ensembl translated-protein identifier for the selected transcript; not the gene symbol. |
| `protein_hgvs` | FGFR, HER2, ERK2 | Protein HGVS notation submitted for annotation, including protein version and three-letter amino-acid substitution. |
| `sift_prediction` | FGFR, HER2, ERK2 | SIFT categorical effect call from Ensembl VEP, e.g. deleterious or tolerated (including low-confidence categories if present). This does not directly predict kinase activation. |
| `sift_score` | FGFR, HER2, ERK2 | Raw SIFT score on 0–1 scale; lower means more likely deleterious. Stored as returned, not an activation probability; filtering uses categorical calls. |
| `polyphen_prediction` | FGFR, HER2, ERK2 | PolyPhen-2 categorical effect call from Ensembl VEP, e.g. benign, possibly_damaging or probably_damaging. Not an activation label. |
| `polyphen_score` | FGFR, HER2, ERK2 | Raw PolyPhen-2 score on 0–1 scale; higher means more likely damaging. Not an activation probability; filtering uses categorical calls. |
| `annotation_status` | FGFR, HER2, ERK2 | OK when sequence/substitution-matched annotation is accepted; ERR-prefixed values denote unavailable, mismatched or conflicting annotation. All three complete tables have OK annotations. |
| `sift_pp_keep` | FGFR, HER2, ERK2 | Boolean annotation-only exclusion-rule retention: annotation_status=OK AND SIFT is not tolerated/tolerated_low_confidence AND PolyPhen is not benign. Scenario eligibility additionally requires scorable=1. |
| `exclusion_reason` | FGFR, HER2, ERK2 | Reason for failing the annotation-only exclusion rule; SIFT_tolerated takes reporting precedence over PolyPhen_benign if both apply. Empty means no exclusion reason; ERR reasons denote annotation problems if present. Does not indicate exclusion from Original scoring. |
| `annotation_origin` | FGFR, HER2, ERK2 | Annotation provenance; these three complete tables use fresh_Ensembl_VEP, not paper-supplied annotations. |
| `retrieved_utc` | FGFR, HER2, ERK2 | UTC retrieval timestamp recorded by the annotation run, not the date of this copy. |
| `Original_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to Original evaluation, otherwise 0. All scorable rows; binary call is raw model probability >= 0.5. |
| `SIFT_PP_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to SIFT_PP evaluation, otherwise 0. Only scorable rows with sift_pp_keep=True; raw binary call. |
| `GT_HS_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to GT_HS evaluation, otherwise 0. All scorable rows; force the call to 1 at GT_HS_flag=1. |
| `SIFT_PP_GT_HS_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to SIFT_PP_GT_HS evaluation, otherwise 0. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. |
| `GBDT_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from gradient-boosted decision-tree model (GBDT), on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `GBDT_prediction` | FGFR, HER2, ERK2 | Baseline binary call for gradient-boosted decision-tree model (GBDT): 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `NGB_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from NGBoost model, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `NGB_prediction` | FGFR, HER2, ERK2 | Baseline binary call for NGBoost model: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `SVM_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from support-vector-machine classifier, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `SVM_prediction` | FGFR, HER2, ERK2 | Baseline binary call for support-vector-machine classifier: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `MNN_member_1_logit` | FGFR, HER2, ERK2 | Raw unbounded pre-sigmoid output of neural checkpoint 1 (of five), not a probability. |
| `MNN_member_2_logit` | FGFR, HER2, ERK2 | Raw unbounded pre-sigmoid output of neural checkpoint 2 (of five), not a probability. |
| `MNN_member_3_logit` | FGFR, HER2, ERK2 | Raw unbounded pre-sigmoid output of neural checkpoint 3 (of five), not a probability. |
| `MNN_member_4_logit` | FGFR, HER2, ERK2 | Raw unbounded pre-sigmoid output of neural checkpoint 4 (of five), not a probability. |
| `MNN_member_5_logit` | FGFR, HER2, ERK2 | Raw unbounded pre-sigmoid output of neural checkpoint 5 (of five), not a probability. |
| `MNN_logit` | FGFR, HER2, ERK2 | Arithmetic mean of the five checkpoint logits, before sigmoid. |
| `MNN_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from five-checkpoint neural-network ensemble, on 0–1 scale. Computed as sigmoid(mean of five member logits), not mean of member probabilities. Never replaced by hotspot/zero overrides. |
| `MNN_prediction` | FGFR, HER2, ERK2 | Baseline binary call for five-checkpoint neural-network ensemble: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `Classical_ensemble_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from arithmetic mean of GBDT, NGB and SVM probabilities, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `GBDT_MNN_ensemble_probability` | FGFR, HER2, ERK2 | Unmodified activating-class probability from arithmetic mean of GBDT and MNN probabilities, on 0–1 scale. Never replaced by hotspot/zero overrides. |
| `Classical_ensemble_prediction` | FGFR, HER2, ERK2 | Baseline binary call for arithmetic mean of GBDT, NGB and SVM probabilities: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `GBDT_MNN_ensemble_prediction` | FGFR, HER2, ERK2 | Baseline binary call for arithmetic mean of GBDT and MNN probabilities: 1 if its raw probability >=0.5, else 0; no annotation or hotspot override. |
| `zero_AND_flag` | FGFR, HER2, ERK2 | 1 if annotation_status=OK AND SIFT is tolerated/tolerated_low_confidence AND PolyPhen is benign; flags a forced-negative call in Zero_AND scenarios, not row removal. |
| `zero_OR_flag` | FGFR, HER2, ERK2 | 1 if annotation_status=OK AND either SIFT is tolerated/tolerated_low_confidence OR PolyPhen is benign; flags a forced-negative call in Zero_OR scenarios, not row removal. |
| `annotation_missing_or_conflicting` | FGFR, HER2, ERK2 | 1 when annotation_status is not OK; otherwise 0. Such rows are retained by Original and full-population zero-override rules but not the exclusion scenario (if scorable). |
| `GBDT_Original_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Original_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Original_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Original_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Original_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Original_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under Original. All scorable rows; binary call is raw model probability >= 0.5. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_SIFT_PP_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under SIFT_PP. Only scorable rows with sift_pp_keep=True; raw binary call. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under GT_HS. All scorable rows; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_SIFT_PP_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under SIFT_PP_GT_HS. Only scorable rows with sift_pp_keep=True; force the call to 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_AND_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to Zero_AND evaluation, otherwise 0. All scorable rows; force the call to 0 at zero_AND_flag=1. |
| `GBDT_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_AND_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under Zero_AND. All scorable rows; force the call to 0 at zero_AND_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_AND_GT_HS_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to Zero_AND_GT_HS evaluation, otherwise 0. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. |
| `GBDT_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_AND_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under Zero_AND_GT_HS. All scorable rows; apply zero_AND_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_OR_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to Zero_OR evaluation, otherwise 0. All scorable rows; force the call to 0 at zero_OR_flag=1. |
| `GBDT_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_OR_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under Zero_OR. All scorable rows; force the call to 0 at zero_OR_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Zero_OR_GT_HS_eligible` | FGFR, HER2, ERK2 | 1 if row contributes to Zero_OR_GT_HS evaluation, otherwise 0. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. |
| `GBDT_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for gradient-boosted decision-tree model (GBDT) under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `NGB_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for NGBoost model under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `SVM_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for support-vector-machine classifier under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `MNN_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for five-checkpoint neural-network ensemble under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `Classical_ensemble_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT, NGB and SVM probabilities under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `GBDT_MNN_ensemble_Zero_OR_GT_HS_prediction` | FGFR, HER2, ERK2 | Binary call for arithmetic mean of GBDT and MNN probabilities under Zero_OR_GT_HS. All scorable rows; apply zero_OR_flag first, then force 1 at GT_HS_flag=1. Empty if ineligible. Experimental labels and raw probabilities are not changed. |
| `activation_conservative` | HER2 | HER2 binary experimental GOF label requiring concordant activity calls in both technical replicates; primary endpoint. |
| `activation_permissive` | HER2 | HER2 broader supplied binary experimental GOF label, including single-replicate calls; equals group_code. |
| `permissive_only_ambiguous` | HER2 | 1 when HER2 is permissive-positive but conservative-negative; otherwise 0. These rows remain negatives in the primary evaluation. |
| `source_ID` | ERK2 | Original supplied spreadsheet identifier, e.g. ERK2_F19W, before normalization to MAPK1. |
| `source_row_1_based` | ERK2 | Original Excel row number including the header (first data row is 2). |
| `fidelity` | ERK2 | Fidelity metadata copied as an integer from the original spreadsheet Fidelity column; not a model confidence measure and not re-derived here. |
| `kd_start` | ERK2 | Inclusive 1-based start of the selected inference KD: ERK2 19. |
| `kd_end` | ERK2 | Inclusive 1-based end of the selected inference KD: ERK2 353. |
| `kd_position_1_based` | ERK2 | Full-protein position minus kd_start plus 1. For unscorable ERK2 rows this can be <=0 or beyond KD length; it is not a valid scored index in those rows. |
| `core_kd_start` | ERK2 | Inclusive 1-based internal database core-KD start: ERK2 25. |
| `core_kd_end` | ERK2 | Inclusive 1-based internal database core-KD end: ERK2 313. |
| `in_internal_core` | ERK2 | 1 if position lies within core_kd_start..core_kd_end, otherwise 0; distinct from membership of the selected extended ERK2 inference KD. |
| `domain_definition` | ERK2 | Selected-domain policy: supplied_extended_kinase_domain for ERK2. |
| `uniprot_accession` | ERK2 | Sequence reference accession: P28482 for MAPK1/ERK2. |

## Data integrity

The populations, schemas and correspondence of each scored file to its complete file
were verified from the CSVs in this directory. This README correction does not change
the data, labels, predictions or analyses. The original spreadsheet `Fidelity` field
is preserved as source metadata and is not reinterpreted as model confidence.

SHA-256 checksums of the six CSVs:

- `ERK2_complete_annotated_predictions.csv`: `dbd70f08888ca6c6d6a09277bfc0a0f8304eaa67c9b13669afcdaa46e3cf4469`
- `ERK2_scored_variants.csv`: `2e425ce141b3becf4082b6fd4177578c1e293dbe23e781494c59e146bc9899d7`
- `FGFR_complete_annotated_predictions.csv`: `24ad7266377d904aa664a8b5cf7cd6ae58a76b5f1f0e99362727b5f42e76dd90`
- `FGFR_scored_variants.csv`: `ca7f8cbbea117afe2fd7f82c92678c3bd3bb99e361bb45118870e49bfa1fbfd4`
- `HER2_complete_annotated_predictions.csv`: `f95febba653459c2dbe847e844745c811e85c4bdea71de19e52bfdc150e2b1c1`
- `HER2_scored_variants.csv`: `19960a97ad923c352e7bbd36fe912fca17fd19b7321cebad3911819ae53e0fc3`
