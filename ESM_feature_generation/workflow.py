"""Validated KD embedding extraction and the original training/PCA data contract.

New extraction is not a claim of bitwise reproduction of historical embeddings.
No function trains a classifier or alters historical input records.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import pickle
import re

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

MODEL_ID = 'facebook/esm2_t33_650M_UR50D'
MODEL_REVISION = '08e4846e537177426273712802403f7ba8261b6c'
FEATURE_DIMS = [1280, 1280, 8, 59, 136, 12]
BOS_COLUMNS = [f'mutant_first_token_ESM_{i}' for i in range(1280)]
RESIDUE_COLUMNS = [f'mutant_ESM_{i}' for i in range(1280)]
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _integer(value, field):
    number = float(value)
    if not np.isfinite(number) or not number.is_integer():
        raise ValueError(f'{field} must be an integer, got {value!r}')
    return int(number)


def validate_variants(variants, references):
    """Preserve source order and labels; mark outside-KD rows unscorable.

    The supplied five identity substitutions remain explicit historical records.
    They are flagged, not silently deleted or relabelled. Validation errors raise.
    """
    required = {'ID', 'Protein', 'Mutation', 'Experimental', 'Fidelity', 'Sequence',
                'mutant_seq', 'KD start', 'KD end', 'KD', 'mutant_KD'}
    if not required.issubset(variants):
        raise ValueError(f'Missing variant columns: {sorted(required - set(variants))}')
    refs = references.rename(columns={'Unnamed: 0': 'Protein'}).copy()
    if not {'Protein', 'Sequence', 'KD start', 'KD end'}.issubset(refs):
        raise ValueError('Reference table lacks protein, sequence or domain coordinates')
    if refs.Protein.isna().any() or refs.Protein.duplicated().any():
        raise ValueError('Reference Protein keys must be unique and nonempty')
    refs = refs.set_index('Protein')
    if variants.ID.isna().any() or variants.ID.duplicated().any():
        raise ValueError('Variant IDs must be unique and nonempty')
    if not variants.Experimental.isin([0, 1]).all():
        raise ValueError('Experimental labels must be binary; no labels are inferred')
    rows = []
    for _, source in variants.iterrows():
        r = source.to_dict()
        gene, mutation = str(r['Protein']), str(r['Mutation'])
        match = re.fullmatch(r'([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])', mutation)
        if match is None or r['ID'] != f'{gene}_{mutation}':
            raise ValueError(f'Invalid or inconsistent substitution/ID: {r["ID"]}')
        if gene not in refs.index:
            raise ValueError(f'No wild-type reference for {gene}')
        reference = refs.loc[gene]
        seq = str(reference.Sequence)
        if not seq or not set(seq).issubset(AMINO_ACIDS):
            raise ValueError(f'Invalid wild-type sequence for {gene}')
        # These are the explicitly selected domains in the supplied project table.
        domain = 'KD2' if gene in {'JAK3', 'RPS6KA5'} else 'KD'
        start = _integer(reference[f'{domain} start'], 'reference KD start')
        end = _integer(reference[f'{domain} end'], 'reference KD end')
        if not 1 <= start <= end <= len(seq):
            raise ValueError(f'Invalid domain bounds for {gene}')
        if (_integer(r['KD start'], 'KD start'), _integer(r['KD end'], 'KD end')) != (start, end):
            raise ValueError(f'Source/reference KD bounds disagree for {r["ID"]}')
        wt, position, mutant = match[1], int(match[2]), match[3]
        if not 1 <= position <= len(seq) or seq[position - 1] != wt:
            raise ValueError(f'WT residue/position mismatch for {r["ID"]}')
        mutated = seq[:position - 1] + mutant + seq[position:]
        if r['Sequence'] != seq or r['mutant_seq'] != mutated:
            raise ValueError(f'Full sequence mismatch for {r["ID"]}')
        if r['KD'] != seq[start - 1:end] or r['mutant_KD'] != mutated[start - 1:end]:
            raise ValueError(f'KD sequence mismatch for {r["ID"]}')
        inside = start <= position <= end
        for field, expected in [('mutation_loc', position), ('kd_mutation_loc', position - start + 1),
                                ('in_KD', int(inside))]:
            if field in r and _integer(r[field], field) != expected:
                raise ValueError(f'{field} disagrees with verified coordinates for {r["ID"]}')
        if 'WT_AA' in r and r['WT_AA'] != wt:
            raise ValueError(f'WT_AA disagrees for {r["ID"]}')
        r.update(position=position, kd_position_1_based=position - start + 1,
                 scorable=int(inside), identity_substitution=int(wt == mutant),
                 inference_status='ready' if inside else 'unscorable_outside_selected_KD')
        rows.append(r)
    return pd.DataFrame(rows)


class ESMExtractor:
    """Pinned HF ESM-2 final layer (33), untruncated KD-only, one sequence at a time."""
    def __init__(self, device=None, local_files_only=False):
        import torch
        from transformers import AutoModel, AutoTokenizer
        self.device = torch.device(device or ('cuda' if torch.cuda.is_available() else 'cpu'))
        self.tokenizer = AutoTokenizer.from_pretrained(
            MODEL_ID, revision=MODEL_REVISION, local_files_only=local_files_only, do_lower_case=False)
        self.model = AutoModel.from_pretrained(
            MODEL_ID, revision=MODEL_REVISION, local_files_only=local_files_only,
            attn_implementation='eager', add_pooling_layer=False).float().eval().to(self.device)

    @property
    def signature(self):
        import torch
        return {'model': MODEL_ID, 'revision': MODEL_REVISION, 'layer': 33,
                'implementation': 'transformers', 'pooling': False, 'dtype': 'float32', 'device': str(self.device),
                'torch': torch.__version__, 'transformers': importlib.metadata.version('transformers')}

    def extract(self, sequence, kd_position, mutant_amino_acid):
        import torch
        if not 1 <= kd_position <= len(sequence) or sequence[kd_position - 1] != mutant_amino_acid:
            raise ValueError('KD mutation index/residue mismatch')
        tokens = self.tokenizer(sequence, return_tensors='pt', truncation=False)
        if tokens['input_ids'].shape[1] != len(sequence) + 2:
            raise ValueError('Unexpected tokenization; expected BOS + residues + EOS')
        token = self.tokenizer.convert_ids_to_tokens(int(tokens['input_ids'][0, kd_position]))
        if token != mutant_amino_acid:
            raise ValueError('Tokenizer residue index does not match the requested mutation')
        tokens = {k: v.to(self.device) for k, v in tokens.items()}
        with torch.inference_mode():
            hidden = self.model(**tokens).last_hidden_state[0]
        # BOS still present: a 1-based residue position equals its token-array index.
        bos = hidden[0].detach().cpu().numpy().astype(np.float32)
        residue = hidden[kd_position].detach().cpu().numpy().astype(np.float32)
        if bos.shape != (1280,) or residue.shape != (1280,) or not np.isfinite([bos, residue]).all():
            raise ValueError('Invalid ESM feature vector')
        return bos, residue


def generate_raw_features(validated, extractor, output_dir):
    """Return all source rows, blank raw vectors for unscorable rows; safe resume."""
    output_dir = Path(output_dir)
    cache = output_dir / 'checkpoints'
    cache.mkdir(parents=True, exist_ok=True)
    signature = extractor.signature
    rows = []
    for _, source in validated.iterrows():
        r = source.to_dict()
        r.update({c: np.nan for c in BOS_COLUMNS + RESIDUE_COLUMNS})
        if r['scorable']:
            seq, index = r['mutant_KD'], int(r['kd_position_1_based'])
            provenance = {'ID': r['ID'], 'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                          'residue_index': index, 'extractor': signature}
            metadata = json.dumps(provenance, sort_keys=True)
            key = hashlib.sha256(metadata.encode()).hexdigest()
            path = cache / f'{key}.npz'
            if path.exists():
                with np.load(path, allow_pickle=False) as saved:
                    if str(saved['metadata']) != metadata:
                        raise ValueError('Embedding cache provenance mismatch')
                    bos, residue = saved['bos'], saved['residue']
            else:
                bos, residue = extractor.extract(seq, index, r['Mutation'][-1])
                temporary = path.with_suffix('.partial.npz')
                np.savez_compressed(temporary, bos=bos, residue=residue, metadata=np.array(metadata))
                temporary.replace(path)
            if bos.shape != (1280,) or residue.shape != (1280,) or not np.isfinite([bos, residue]).all():
                raise ValueError(f'Invalid cached vectors for {r["ID"]}')
            r.update(zip(BOS_COLUMNS, bos.tolist()))
            r.update(zip(RESIDUE_COLUMNS, residue.tolist()))
            r['inference_status'] = 'embedded'
        rows.append(r)
    result = pd.DataFrame(rows)
    if result.ID.tolist() != validated.ID.tolist():
        raise ValueError('Output coverage/order changed')
    result.to_csv(output_dir / 'KD_raw_mutant_ESM_features.csv', index=False)
    (output_dir / 'embedding_provenance.json').write_text(json.dumps({
        'extractor': signature, 'rows': len(result), 'embedded': int(result.scorable.sum()),
        'historical_training_features_replaced': False}, indent=2))
    return result


def load_training_contract(csv_path, pickle_path):
    """Read only trusted project pickles: pickle can execute code when loaded."""
    # Preserve the decimal-to-float round trip for exact CSV/pickle comparisons.
    table = pd.read_csv(csv_path, float_precision='round_trip')
    with Path(pickle_path).open('rb') as stream:
        package = pickle.load(stream)
    if len(package) != 7:
        raise ValueError('Expected seven processed-data items')
    features, labels, fidelity, privileged, first_basis, residue_basis, dims = package
    if list(dims) != FEATURE_DIMS or features.shape != (len(table), sum(FEATURE_DIMS)):
        raise ValueError('Unexpected original training dimensions')
    if len(table) != 1003 or table.shape[1] != 2779:
        raise ValueError('This validator expects the unchanged 1003-row original training dataset')
    if table.columns[1:2561].tolist() != BOS_COLUMNS + RESIDUE_COLUMNS:
        raise ValueError('Raw ESM column layout differs from original contract')
    if table.columns[2561:2569].tolist() != [f'Kinase_group_one_hot{i}' for i in range(8)]:
        raise ValueError('Kinase-family feature order differs')
    if table.columns[-3:].tolist() != ['Label', 'Kinase', 'mutation']:
        raise ValueError('Expected separate label and identifier columns')
    if not np.allclose(table.iloc[:, 1:2776].to_numpy(float), features, rtol=0, atol=0, equal_nan=True):
        raise ValueError('CSV and processed feature values/order differ')
    if not np.array_equal(table.Label.to_numpy(), labels):
        raise ValueError('CSV and processed labels differ')
    if not np.isin(labels, [0, 1]).all() or int(np.sum(labels)) != 329:
        raise ValueError('Original labels differ')
    if not all(np.asarray(a).shape == (len(table),) for a in [labels, fidelity, privileged]):
        raise ValueError('Metadata vectors must align to every training row')
    if not all(np.isin(a, [0, 1]).all() for a in [fidelity, privileged]):
        raise ValueError('Expected binary fidelity and PI vectors')
    if first_basis.shape != (196, 1280) or residue_basis.shape != (53664, 1280):
        raise ValueError('Reference bases differ from the authoritative original artifact')
    if not np.isfinite(features[:, :2560]).all():
        raise ValueError('Historical ESM blocks must be complete; do not drop the 17 source out-of-KD rows')
    return table, package


def fit_reference_pcas(first_basis, residue_basis):
    """Match original trainers: unscaled reference embeddings; NOT fitting on new cases."""
    for basis in [first_basis, residue_basis]:
        if basis.ndim != 2 or basis.shape[1] != 1280 or len(basis) < 50 or not np.isfinite(basis).all():
            raise ValueError('PCA basis must be finite, N x 1280, with N >= 50')
    first = PCA(n_components=50, random_state=18).fit(first_basis)
    residue = PCA(n_components=50, random_state=18).fit(residue_basis)
    return first, residue


def assemble_training_features(features, first_pca, residue_pca):
    if features.ndim != 2 or features.shape[1] != sum(FEATURE_DIMS):
        raise ValueError('Expected original 2775-column raw feature matrix')
    if first_pca.components_.shape != (50, 1280) or residue_pca.components_.shape != (50, 1280):
        raise ValueError('Expected saved 50-component PCAs for each raw embedding block')
    return np.hstack([first_pca.transform(features[:, :1280]),
                      residue_pca.transform(features[:, 1280:2560]), features[:, 2560:]])


def apply_saved_classical_preprocessing(features, first_pca, residue_pca, imputer=None, scaler=None):
    """Inference only: transform with existing artifacts, never fit on external cases."""
    selected = assemble_training_features(features, first_pca, residue_pca)[:, :167]
    if imputer is not None:
        selected = imputer.transform(selected)
    if scaler is not None:
        selected = scaler.transform(selected)
    return selected


def map_original_training_ids(training_table, variants):
    """Explicit legacy-to-source mapping, not a change to stored historical IDs."""
    keys = training_table[['Kinase', 'mutation', 'Label']].copy()
    keys['training_row'] = np.arange(len(keys))
    keys['Protein'] = keys.Kinase.replace({'BRAF-NEW': 'BRAF', 'BRAF-OLD': 'BRAF',
                                          'MEK1': 'MAP2K1', 'ERBB2': 'HER2'})
    keys['Mutation'] = keys.mutation
    egfr = keys.Kinase.eq('EGFR')
    keys.loc[egfr, 'Mutation'] = keys.loc[egfr, 'mutation'].str.replace(
        r'(\d+)', lambda match: str(int(match[1]) + 24), regex=True)
    mapped = keys.merge(variants[['ID', 'Protein', 'Mutation', 'Experimental']],
                        on=['Protein', 'Mutation'], how='left', validate='one_to_one', sort=False)
    if mapped.ID.isna().any() or not mapped.ID.is_unique or not mapped.Label.eq(mapped.Experimental).all():
        raise ValueError('Historical/source identifier mapping or labels disagree')
    if mapped.training_row.tolist() != list(range(len(training_table))):
        raise ValueError('Historical row order changed')
    return mapped
