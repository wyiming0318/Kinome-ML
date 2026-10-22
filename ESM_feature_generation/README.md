# ESM feature generation and the training/PCA workflow

This folder contains **code and original sequence inputs**, not a replacement
training dataset. The previous exploratory notebook and its four embedding/PCA
exports were removed because they did not faithfully implement the project's
KD feature and training-data contracts. Their removal does not establish that
the original models were trained on those exports; they remain in Git history.
The license and both supplied source CSVs are preserved unchanged.

## Run the notebooks

Use Python 3.11 in a new environment, install `requirements.txt`, and register
that environment as a Jupyter kernel. Open notebooks with this folder as the
working directory. The requirements record the tested notebook environment,
**not** the historical classifier-training environment. Loading old serialized
classifiers may require their original dependency versions.

1. **`01_KD_ESM_Feature_Generation.ipynb`** validates `mutant_df.csv` against
   `wildtype_kinase_data.csv`, then optionally generates KD-only raw ESM features.
   It does not run a classifier or fit PCA. Embedding generation is deliberately
   opt-in because the 650M-parameter model and complete extraction are expensive.
2. **`02_Training_Data_and_PCA_Workflow.ipynb`** validates externally supplied
   original training data and demonstrates the exact feature-block ordering and
   PCA fitting configuration used in the trainers. Full reference-PCA fitting
   is also opt-in. Use saved model-associated PCA objects for existing inference.

`workflow.py` implements the shared validation/transformation logic;
`python -m pytest -q test_workflow.py` runs offline tests without downloading
the model. No classifier is trained by these notebooks or tests.

## Inputs and optional execution settings

| Input / setting | Meaning |
|---|---|
| `mutant_df.csv` | Original 1,003 source records with IDs, experimental labels, fidelity, reference/mutated full sequences and KD slices; not the full model feature matrix. |
| `wildtype_kinase_data.csv` | 196 wild-type reference sequences and domain coordinates. JAK3 and RPS6KA5 explicitly use the second domain. |
| `KINOME_SOURCE_DIR` | Optional directory containing the two source CSVs; defaults to this folder. |
| `KINOME_OUTPUT_DIR` | Local generated-output directory; defaults to ignored `outputs/`. Inputs are never overwritten. |
| `KINOME_RUN_ESM=1` | Enable actual ESM extraction; unset by default. First use can download the pinned HF model. |
| `KINOME_ESM_LIMIT=N` | Optional first-N-row smoke run, retaining input order. Default 0 means all rows; a limited export is not a complete dataset. |
| `KINOME_LOCAL_MODEL_ONLY=1` | Require the model already to be cached; no model download. |
| `KINOME_TRAINING_CSV` | Original full training CSV, default `../data/data.csv`; the investigator's local equivalent is `dataset.csv`. Not `mutant_df.csv`. |
| `KINOME_DATASET_PKL` | Original seven-item processed-data pickle, default `../processed_data/dataset.pkl`. |
| `KINOME_TRUST_PICKLE=1` | Explicit confirmation that the supplied pickle is a trusted project file; pickle loading can execute code. |
| `KINOME_RUN_PCA=1` | Fit the two reference PCAs for demonstration after validation; unset by default. Never use this to refit preprocessing for external inference. |

Supply the original training files separately using the main README's data
download or your project copy. The notebooks do not reconstruct missing
historical features from the student's CSV exports or from new embeddings.

## Sequence and token contract

Full-protein mutation positions and KD bounds are **1-based, inclusive**.
Each source WT residue, mutated sequence and KD slice is checked. Embeddings
use only `mutant_KD`, never the full protein: ESM-2
`facebook/esm2_t33_650M_UR50D`, final layer 33, width 1,280, float32. The HF
revision is pinned in `workflow.py`; CPU/CUDA selection is automatic.

Token zero is BOS. While BOS is present, the residue token index equals
`full_position - kd_start + 1`. If BOS/EOS are stripped, subtract one from
that KD position. The code also verifies the token's amino-acid identity.
Only the mutant BOS and mutant mutation-residue vectors are sample features;
WT vectors and WT-minus-mutant differences are not appended to this contract.

All 1,003 source rows and labels are retained. The 17 mutations outside the
source KD definitions are unscorable for **new KD-only residue extraction**;
their generated feature cells remain missing. Five supplied identity
substitutions are flagged as historical no-op records, not silently removed;
in-KD identity records may be embedded but are not genuine missense changes.
In contrast, the existing historical training matrix contains complete raw
features for all 1,003 records and is preserved as supplied. Do not overwrite
it with, or inner-join it down to, the 986 in-KD source rows.

New extraction writes `outputs/raw_esm/KD_raw_mutant_ESM_features.csv`, with
the original source columns, validation/status metadata, and
`mutant_first_token_ESM_0..1279` / `mutant_ESM_0..1279` raw columns. These names
match the raw training blocks, but this export is **not** a complete training
table: the 67 basic/family features and structural features are not fabricated.
Per-record caches include sequence, token index, model and runtime provenance.
Errors fail explicitly; reruns only reuse compatible completed caches.

## Original training layout and PCA

The original CSV has 1,003 rows and 2,779 columns, including one saved dataframe
index, 2,775 features, and `Label`, `Kinase`, `mutation`. It has 329 positives
and 674 negatives. The seven-item processed pickle is:

`[features, labels, fidelity, PI_flags, BOS_reference_basis, residue_reference_basis, feature_dimensions]`.

| Feature block in stored order | Raw width | After ESM PCA |
|---|---:|---:|
| Mutant BOS ESM | 1,280 | 50 |
| Mutant residue ESM | 1,280 | 50 |
| Kinase-family one-hot | 8 | 8 |
| Basic sequence/physicochemical | 59 | 59 |
| MD | 136 | 136 |
| Docking | 12 | 12 |
| Total feature matrix | 2,775 | 315 |

Two `PCA(n_components=50, random_state=18)` objects are fitted directly on the
supplied **196 x 1,280 BOS reference matrix** and **53,664 x 1,280 residue
reference matrix**, respectively. No `StandardScaler` precedes PCA. These
bases are separate from the sample matrix: they are neither the 1,003 mutant
sample vectors nor an automatically regenerated WT table. They must be
supplied intact. The original trainers fit reference PCA before cross-validation,
not separately in each fold. The 100-dimensional combined ESM input means
**50 + 50**, not 100 per embedding. New PCA results can differ across library
versions; fitting this demonstration does not replace frozen model transforms.

The 8 family indicators are not the experimental label. `PI_flags` mark the
215 records with privileged structural features (788 without); fidelity is a
separate binary vector. The shown original trainer paths load fidelity but do
not use it as a sample weight. Missing structural values remain missing.

Historical identifiers are not simply interchangeable with source identifiers:
the explicit mapping covers BRAF-NEW/BRAF-OLD to BRAF, MEK1 to MAP2K1,
ERBB2 to HER2, and the original EGFR +24 source-coordinate offset. That offset
is specific to this original training/source mapping, not a rule for new EGFR
variants. Existing IDs, source labels and training row order remain unchanged.

## Model-specific use (not classifier training)

- **Classical no-PI models:** 50 BOS PCs + 50 residue PCs + 8 family + 59 basic
  features = 167. Apply model-associated imputation (if enabled), then scaling,
  then classifier. In CV, the imputer/scaler are fitted to the training fold;
  inference only calls their saved `transform` methods.
- **Original `NNtrainer.py`:** its current `use_basic=False` branch actually
  selects a 167-wide main input; do not interpret the flag name literally.
  Its privileged/joint branch includes 148 structural features plus 100 ESM
  PCs (248 columns). The unscaled trainer does not standardize those inputs.
- **`NNtrainerScaled.py` local variant:** same PCA configuration; scaling is
  after PCA/feature assembly, not before raw-embedding PCA.
- **Released collaborator checkpoints used in the current inference GUI:**
  `sub_shared_dim=100`, so the actual main input is only the two ESM PCA blocks.
  Read checkpoint dimensions rather than assuming the original trainer flags
  describe every released artifact. This notebook does not publish weights.

Existing inference must reuse **each model's own** PCA, imputer and scaler;
do not fit on a new test dataset or assume all released PCA artifacts are
interchangeable. Original trainers use repeated five-fold stratified CV (five
repeats); fitted final-model training accuracy is not an OOF estimate.

## Provenance boundary

The supplied historical training features and frozen transforms are authoritative
for reproducing existing results. This corrected extraction implementation
matches the intended KD/token contract, **not proven bitwise historical
embedding provenance**. The superseded exports did not match the original
training vectors; their KD first-token mutant vectors duplicated the whole-
protein vectors. The original extraction implementation and exact reference-
basis construction are not inferred from those defective files. No historical
training record, label, published metric or saved classifier is changed here.
