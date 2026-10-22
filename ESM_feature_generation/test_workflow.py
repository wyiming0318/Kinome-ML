"""Run: python -m pytest -q test_workflow.py (no network/model downloads)."""
from pathlib import Path
from types import SimpleNamespace
import tempfile

import numpy as np
import pandas as pd
import pytest
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

from workflow import (BOS_COLUMNS, RESIDUE_COLUMNS, ESMExtractor, validate_variants,
                      generate_raw_features, fit_reference_pcas, assemble_training_features,
                      apply_saved_classical_preprocessing, map_original_training_ids)


def fixture(mutations=('D3A', 'F5Y', 'C2A', 'G6W', 'E4E')):
    sequence = 'ACDEFGH'
    rows = []
    for mutation in mutations:
        pos = int(mutation[1:-1])
        mutated = sequence[:pos-1] + mutation[-1] + sequence[pos:]
        rows.append(dict(ID='TEST_'+mutation, Protein='TEST', Mutation=mutation,
                         Experimental=1, Fidelity=0, Sequence=sequence, mutant_seq=mutated,
                         KD=sequence[2:5], mutant_KD=mutated[2:5], **{'KD start': 3, 'KD end': 5}))
    refs = pd.DataFrame([{'Protein': 'TEST', 'Sequence': sequence, 'KD start': 3, 'KD end': 5}])
    return pd.DataFrame(rows), refs


def test_inclusive_bounds_order_and_identity():
    variants, refs = fixture()
    result = validate_variants(variants, refs)
    assert result.ID.tolist() == variants.ID.tolist()
    assert result.scorable.tolist() == [1, 1, 0, 0, 1]
    assert result.kd_position_1_based.tolist() == [1, 3, 0, 4, 2]
    assert result.identity_substitution.tolist() == [0, 0, 0, 0, 1]
    pd.testing.assert_series_equal(result.Experimental, variants.Experimental)


@pytest.mark.parametrize('field,value', [('Mutation', 'A3V'), ('Sequence', 'XXXXXXX'),
    ('mutant_seq', 'ACDEFGH'), ('mutant_KD', 'XXX'), ('KD end', 4),
    ('Experimental', 2), ('ID', 'wrong'), ('Mutation', 'D3*'), ('Mutation', 'D30A')])
def test_invalid_records_fail(field, value):
    variants, refs = fixture(('D3A',))
    variants.loc[0, field] = value
    with pytest.raises(ValueError):
        validate_variants(variants, refs)


def test_duplicate_ids_and_unknown_reference():
    variants, refs = fixture(('D3A', 'D3A'))
    with pytest.raises(ValueError):
        validate_variants(variants, refs)
    variants, refs = fixture(('D3A',))
    refs.loc[0, 'Protein'] = 'OTHER'
    with pytest.raises(ValueError):
        validate_variants(variants, refs)


@pytest.mark.parametrize('gene', ['JAK3', 'RPS6KA5'])
def test_second_domain_selected(gene):
    variants, refs = fixture(('D3A',))
    variants['Protein'] = gene
    variants['ID'] = gene + '_D3A'
    refs['Protein'] = gene
    refs['KD2 start'], refs['KD2 end'] = 3, 5
    refs['KD start'], refs['KD end'] = 1, 2
    assert validate_variants(variants, refs).scorable.iloc[0] == 1


def test_token_index_with_bos_present():
    import torch
    class Tokenizer:
        def __call__(self, seq, **kwargs):
            self.seq = seq
            return {'input_ids': torch.tensor([[0]+list(range(1,len(seq)+1))+[0]])}
        def convert_ids_to_tokens(self, idx):
            return self.seq[idx-1]
    class Model:
        def __call__(self, **kwargs):
            count = kwargs['input_ids'].shape[1]
            return SimpleNamespace(last_hidden_state=torch.arange(count).reshape(1,count,1).expand(1,count,1280).float())
    extractor = ESMExtractor.__new__(ESMExtractor)
    extractor.device, extractor.tokenizer, extractor.model = torch.device('cpu'), Tokenizer(), Model()
    for position in [1, 3]:
        bos, residue = extractor.extract('DEF', position, 'DEF'[position-1])
        assert (bos == 0).all() and (residue == position).all()
    with pytest.raises(ValueError):
        extractor.extract('DEF', 2, 'A')


def test_embedding_coverage_cache_and_no_silent_failure(tmp_path):
    class Extractor:
        signature = {'model': 'unit-test', 'layer': 33}
        calls = 0
        def extract(self, sequence, index, aa):
            self.calls += 1
            return np.full(1280, len(sequence)), np.full(1280, index)
    variants, refs = fixture()
    verified = validate_variants(variants, refs)
    extractor = Extractor()
    result = generate_raw_features(verified, extractor, tmp_path)
    assert extractor.calls == 3
    assert result.loc[result.scorable.eq(0), BOS_COLUMNS+RESIDUE_COLUMNS].isna().all().all()
    assert result.loc[0, 'mutant_ESM_0'] == 1 and result.loc[1, 'mutant_ESM_0'] == 3
    generate_raw_features(verified, extractor, tmp_path)
    assert extractor.calls == 3
    # Changing provenance must not reuse an existing cached embedding.
    extractor.signature = {'model': 'unit-test-changed', 'layer': 33}
    generate_raw_features(verified, extractor, tmp_path)
    assert extractor.calls == 6
    class Broken(Extractor):
        def extract(self, *args):
            raise RuntimeError('inference failure')
    with pytest.raises(RuntimeError):
        generate_raw_features(verified, Broken(), tmp_path / 'failed')


def test_pca_exact_contract_and_saved_preprocessing():
    rng = np.random.default_rng(18)
    first_basis = rng.normal(size=(64, 1280)) * np.linspace(1, 5, 1280)
    residue_basis = rng.normal(size=(70, 1280))
    first, residue = fit_reference_pcas(first_basis, residue_basis)
    assert first.components_.shape == residue.components_.shape == (50, 1280)
    np.testing.assert_allclose(first.mean_, first_basis.mean(axis=0))
    direct = PCA(n_components=50, random_state=18).fit(first_basis)
    np.testing.assert_allclose(first.components_, direct.components_)
    features = rng.normal(size=(8, 2775))
    features[0, -1] = np.nan
    transformed = assemble_training_features(features, first, residue)
    assert transformed.shape == (8, 315)
    np.testing.assert_allclose(transformed[:,100:], features[:,2560:], equal_nan=True)
    assert transformed[:,:167].shape == (8,167) and transformed[:,:100].shape == (8,100)
    imputer = SimpleImputer().fit(transformed[:,:167])
    scaler = StandardScaler().fit(imputer.transform(transformed[:,:167]))
    actual = apply_saved_classical_preprocessing(features,first,residue,imputer,scaler)
    expected = scaler.transform(imputer.transform(transformed[:,:167]))
    np.testing.assert_allclose(actual,expected)
    with pytest.raises(ValueError):
        fit_reference_pcas(first_basis[:20], residue_basis)


def test_legacy_aliases_and_egfr_offset():
    train = pd.DataFrame({'Kinase':['BRAF-NEW','MEK1','EGFR'], 'mutation':['V600E','P124S','L858R'], 'Label':[1,0,1]})
    variants = pd.DataFrame({'Protein':['BRAF','MAP2K1','EGFR'], 'Mutation':['V600E','P124S','L882R'],
                             'Experimental':[1,0,1], 'ID':['BRAF_V600E','MAP2K1_P124S','EGFR_L882R']})
    mapped = map_original_training_ids(train,variants)
    assert mapped.ID.tolist() == variants.ID.tolist()
    variants.loc[0,'Experimental'] = 0
    with pytest.raises(ValueError):
        map_original_training_ids(train,variants)
