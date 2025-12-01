#!/usr/bin/env python3
"""
Tests for DeepCoy decoy generation pipeline.

Run with: pytest tests/test_deepcoy_generation.py -v
"""

import os
import sys
import json
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import pandas as pd
import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from decoy.deepcoy_generator import (
    DeepCoyConfig,
    MolecularProperties,
    DecoyPropertyMatcher,
    ActiveMoleculeLoader,
    DecoySelector,
    DecoyGenerationPipeline
)
from decoy.decoy_evaluation import (
    auc_one_positive,
    calculate_doe_score,
    calculate_lads_score,
    evaluate_decoys
)


# Test data
SAMPLE_ACTIVES = [
    "O=C1c2ccccc2C(=O)c2c1cc(S(=O)(=O)[O-])c(O)c2O.[Na+]",  # Anthraquinone derivative
    "CC(C)=CC[C@@H](O)C1=CC(=O)c2c(O)ccc(O)c2C1=O",  # Shikonin
    "O=C1C=C(CSC(=S)N2CCN(c3ccccc3)CC2)C(=O)c2ccccc21",  # Naphthoquinone
]

SAMPLE_DECOYS = [
    "c1ccc2c(c1)ccc1ccccc12",  # Naphthalene
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CN1CCN(c2ccc(C(=O)N3CCN(c4ccccn4)CC3)cc2)CC1",  # Drug-like
]


class TestMolecularProperties:
    """Tests for MolecularProperties class."""
    
    def test_from_smiles_valid(self):
        """Test property calculation from valid SMILES."""
        props = MolecularProperties.from_smiles("CCO")  # Ethanol
        
        assert props is not None
        assert props.smiles == "CCO"
        assert 40 < props.mol_weight < 50  # ~46.07
        assert -1 < props.logp < 0.5
        assert props.hba >= 0
        assert props.hbd >= 0
    
    def test_from_smiles_invalid(self):
        """Test handling of invalid SMILES."""
        props = MolecularProperties.from_smiles("invalid_smiles")
        assert props is None
    
    def test_to_dict(self):
        """Test conversion to dictionary."""
        props = MolecularProperties.from_smiles("CCO")
        d = props.to_dict()
        
        assert 'smiles' in d
        assert 'MW' in d
        assert 'logP' in d
        assert 'HBA' in d
        assert 'HBD' in d


class TestDecoyPropertyMatcher:
    """Tests for DecoyPropertyMatcher class."""
    
    @pytest.fixture
    def matcher(self):
        config = DeepCoyConfig()
        return DecoyPropertyMatcher(config)
    
    def test_property_match_score_identical(self, matcher):
        """Test that identical molecules have perfect match score."""
        props = MolecularProperties.from_smiles("CCO")
        score = matcher.property_match_score(props, props)
        assert score == 1.0
    
    def test_property_match_score_different(self, matcher):
        """Test property matching between different molecules."""
        props1 = MolecularProperties.from_smiles("CCO")  # Ethanol
        props2 = MolecularProperties.from_smiles("CCCCCCCCCCCCCCCC")  # Hexadecane
        
        score = matcher.property_match_score(props1, props2)
        assert 0 <= score <= 1
        assert score <= 0.5  # Should be quite different
    
    def test_tanimoto_similarity_identical(self, matcher):
        """Test that identical molecules have similarity 1.0."""
        smiles = "CCO"
        sim = matcher.calculate_tanimoto_similarity(smiles, smiles)
        assert sim == 1.0
    
    def test_tanimoto_similarity_different(self, matcher):
        """Test similarity between different molecules."""
        sim = matcher.calculate_tanimoto_similarity(
            "CCO",  # Ethanol
            "c1ccccc1"  # Benzene
        )
        assert 0 <= sim <= 1
        assert sim < 0.5  # Should be quite different
    
    def test_max_similarity_to_actives(self, matcher):
        """Test finding max similarity to a list of actives."""
        max_sim = matcher.max_similarity_to_actives(
            "CCO",
            ["CCCO", "CCCCO", "c1ccccc1"]
        )
        
        assert 0 <= max_sim <= 1
        # Should be most similar to CCCO (propanol)


class TestActiveMoleculeLoader:
    """Tests for ActiveMoleculeLoader class."""
    
    @pytest.fixture
    def loader(self):
        config = DeepCoyConfig()
        return ActiveMoleculeLoader(config)
    
    def test_load_from_chembl_csv(self, loader):
        """Test loading from ChEMBL CSV format."""
        # Create temporary CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("molecule_chembl_id,canonical_smiles,standard_value\n")
            f.write("CHEMBL1,CCO,100\n")
            f.write("CHEMBL2,CCCO,200\n")
            temp_path = f.name
        
        try:
            df = loader.load_from_chembl_csv(temp_path)
            assert len(df) == 2
            assert 'canonical_smiles' in df.columns
        finally:
            os.unlink(temp_path)
    
    def test_load_from_smiles_file(self, loader):
        """Test loading from SMILES file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.smi', delete=False) as f:
            f.write("CCO\n")
            f.write("CCCO\n")
            f.write("CCCCO\n")
            temp_path = f.name
        
        try:
            df = loader.load_from_smiles_file(temp_path)
            assert len(df) == 3
        finally:
            os.unlink(temp_path)
    
    def test_calculate_properties(self, loader):
        """Test property calculation for dataframe."""
        df = pd.DataFrame({'canonical_smiles': ['CCO', 'CCCO']})
        result = loader.calculate_properties(df)
        
        assert 'MW' in result.columns
        assert 'logP' in result.columns
        assert len(result) == 2


class TestDecoyEvaluation:
    """Tests for decoy evaluation functions."""
    
    def test_auc_one_positive_all_greater(self):
        """Test AUC when positive is greater than all negatives."""
        auc = auc_one_positive(10.0, np.array([1, 2, 3, 4, 5]))
        assert auc == 1.0
    
    def test_auc_one_positive_all_less(self):
        """Test AUC when positive is less than all negatives."""
        auc = auc_one_positive(0.0, np.array([1, 2, 3, 4, 5]))
        assert auc == 0.0
    
    def test_auc_one_positive_ties(self):
        """Test AUC with ties."""
        auc = auc_one_positive(5.0, np.array([5, 5, 5, 5, 5]))
        assert auc == 0.5  # All ties
    
    def test_calculate_doe_score(self):
        """Test DOE score calculation."""
        # Create synthetic property matrices
        active_props = np.array([
            [300, 2.5, 3, 2, 4],  # MW, logP, HBA, HBD, RotB
            [320, 2.8, 2, 1, 3],
        ])
        
        decoy_props = np.array([
            [310, 2.6, 3, 2, 4],  # Similar properties
            [315, 2.7, 2, 2, 3],
            [305, 2.4, 3, 1, 4],
        ])
        
        doe, per_prop = calculate_doe_score(
            active_props, decoy_props,
            ['MW', 'logP', 'HBA', 'HBD', 'RotB']
        )
        
        assert 0 <= doe <= 0.5
        assert len(per_prop) == 5
    
    def test_calculate_lads_score(self):
        """Test LADS score calculation."""
        actives = ["CCO", "CCCO"]
        decoys = ["c1ccccc1", "c1ccc2ccccc2c1"]  # Very different
        
        lads = calculate_lads_score(actives, decoys)
        
        assert 0 <= lads <= 1
        # Should be low since aromatics are different from alcohols


class TestDecoySelector:
    """Tests for DecoySelector class."""
    
    @pytest.fixture
    def selector(self):
        config = DeepCoyConfig(
            num_decoys_per_active=5,
            max_tanimoto_similarity=0.35,
            min_property_match_score=0.5
        )
        return DecoySelector(config)
    
    def test_select_decoys_filtering(self, selector):
        """Test that selector properly filters candidates."""
        active = "c1ccccc1"  # Benzene
        candidates = [
            "c1ccc2ccccc2c1",  # Naphthalene (similar)
            "CCCCCCCCCCCC",  # Dodecane (very different)
            "c1ccccc1C",  # Toluene (too similar)
        ]
        all_actives = [active]
        
        selected = selector.select_decoys(active, candidates, all_actives)
        
        # Should filter out too similar candidates
        assert len(selected) <= len(candidates)


class TestIntegration:
    """Integration tests for the full pipeline."""
    
    def test_pipeline_process_target(self):
        """Test processing a target with mock data."""
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("molecule_chembl_id,canonical_smiles,standard_value\n")
            for i, smi in enumerate(SAMPLE_ACTIVES[:2]):
                f.write(f"CHEMBL{i},{smi},100\n")
            input_path = f.name
        
        with tempfile.TemporaryDirectory() as output_dir:
            try:
                config = DeepCoyConfig(
                    num_decoys_per_active=5,
                    num_candidates_per_active=10
                )
                pipeline = DecoyGenerationPipeline(config)
                
                result = pipeline.process_target(
                    input_path,
                    output_dir,
                    "test_target"
                )
                
                assert 'num_actives' in result
                assert result['num_actives'] > 0
                
                # Check output files were created
                assert os.path.exists(os.path.join(output_dir, "test_target_actives.csv"))
                assert os.path.exists(os.path.join(output_dir, "test_target_summary.json"))
                
            finally:
                os.unlink(input_path)
    
    def test_full_evaluation(self):
        """Test full evaluation with sample data."""
        # Skip if RDKit not available
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not available")
        
        result = evaluate_decoys(
            SAMPLE_ACTIVES,
            SAMPLE_DECOYS,
            calculate_ml=False  # Skip ML for faster test
        )
        
        assert result.n_actives == len(SAMPLE_ACTIVES)
        assert result.n_decoys == len(SAMPLE_DECOYS)
        assert 0 <= result.doe_score <= 0.5
        assert 0 <= result.lads_score <= 1


class TestDeepCoyConfig:
    """Tests for DeepCoyConfig dataclass."""
    
    def test_default_values(self):
        """Test default configuration values."""
        config = DeepCoyConfig()
        
        assert config.num_decoys_per_active == 50
        assert config.max_tanimoto_similarity == 0.35
        assert config.dataset == "zinc"
    
    def test_custom_values(self):
        """Test custom configuration values."""
        config = DeepCoyConfig(
            num_decoys_per_active=100,
            max_tanimoto_similarity=0.4,
            mw_tolerance=50.0
        )
        
        assert config.num_decoys_per_active == 100
        assert config.max_tanimoto_similarity == 0.4
        assert config.mw_tolerance == 50.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
