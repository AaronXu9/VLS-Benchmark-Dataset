"""
Unit tests for PLI Similarity Scorer.

Tests cover:
1. Core scoring functions (Tanimoto, SuCOS, feature maps)
2. Molecule transformations
3. Foldseek integration
4. End-to-end scoring workflow
"""

import os
import tempfile
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

# Add scripts directory to path to allow importing the module
sys.path.append(str(Path(__file__).parent.parent / "scripts"))

# Import the module under test
from pli_similarity_scorer import (
    PLISimilarityScorer,
    FoldseekRunner,
    FoldseekAlignment,
    SystemInfo,
    SimilarityScore,
    get_feature_map_score,
    get_sucos_score,
    align_molecules,
    calculate_tanimoto_similarity,
    apply_transformation_to_molecule,
)


class TestTanimotoSimilarity(unittest.TestCase):
    """Test 2D Tanimoto fingerprint similarity."""
    
    def test_identical_molecules(self):
        """Identical molecules should have Tanimoto = 1.0."""
        smiles = "CCO"  # Ethanol
        similarity = calculate_tanimoto_similarity(smiles, smiles)
        self.assertAlmostEqual(similarity, 1.0, places=5)
    
    def test_similar_molecules(self):
        """Similar molecules should have high Tanimoto."""
        # Ethanol vs Methanol - similar small alcohols
        similarity = calculate_tanimoto_similarity("CCO", "CO")
        self.assertGreater(similarity, 0.3)
        self.assertLess(similarity, 1.0)
    
    def test_different_molecules(self):
        """Very different molecules should have low Tanimoto."""
        # Ethanol vs benzene - structurally different
        similarity = calculate_tanimoto_similarity("CCO", "c1ccccc1")
        self.assertLess(similarity, 0.5)
    
    def test_invalid_smiles(self):
        """Invalid SMILES should return NaN."""
        similarity = calculate_tanimoto_similarity("invalid", "CCO")
        self.assertTrue(np.isnan(similarity))
    
    def test_drug_molecules(self):
        """Test with drug-like molecules."""
        # Aspirin vs Ibuprofen
        aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        similarity = calculate_tanimoto_similarity(aspirin, ibuprofen)
        self.assertGreater(similarity, 0.0)
        self.assertLess(similarity, 1.0)


class TestFeatureMapScore(unittest.TestCase):
    """Test pharmacophore feature overlap scoring."""
    
    def setUp(self):
        """Create test molecules with 3D conformers."""
        # Simple molecules with pharmacophore features
        self.mol_1 = Chem.MolFromSmiles("CCO")  # Ethanol - has donor/acceptor
        self.mol_2 = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        
        # Add 3D conformers
        AllChem.EmbedMolecule(self.mol_1, randomSeed=42)
        AllChem.EmbedMolecule(self.mol_2, randomSeed=42)
    
    def test_feature_map_score_range(self):
        """Feature map score should be between 0 and 1."""
        score = get_feature_map_score(self.mol_1, self.mol_2)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)
    
    def test_identical_molecules_feature_score(self):
        """Identical molecules should have feature score = 1."""
        mol_copy = Chem.Mol(self.mol_1)
        score = get_feature_map_score(self.mol_1, mol_copy)
        self.assertAlmostEqual(score, 1.0, places=3)


class TestSuCOSScore(unittest.TestCase):
    """Test SuCOS (Shape and Color Overlap Score) calculation."""
    
    def setUp(self):
        """Create test molecules."""
        self.mol_1 = Chem.MolFromSmiles("c1ccccc1O")  # Phenol
        self.mol_2 = Chem.MolFromSmiles("c1ccccc1N")  # Aniline
        
        AllChem.EmbedMolecule(self.mol_1, randomSeed=42)
        AllChem.EmbedMolecule(self.mol_2, randomSeed=42)
    
    def test_sucos_score_range(self):
        """SuCOS score should be between 0 and 1."""
        score = get_sucos_score(self.mol_1, self.mol_2)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)
    
    def test_identical_molecules_sucos(self):
        """Identical molecules should have high SuCOS."""
        mol_copy = Chem.Mol(self.mol_1)
        score = get_sucos_score(self.mol_1, mol_copy)
        self.assertGreater(score, 0.8)
    
    def test_similar_molecules_sucos(self):
        """Similar molecules (phenol vs aniline) should have moderate SuCOS."""
        score = get_sucos_score(self.mol_1, self.mol_2)
        self.assertGreater(score, 0.3)


class TestMoleculeTransformation(unittest.TestCase):
    """Test coordinate transformation functions."""
    
    def setUp(self):
        """Create a test molecule with known coordinates."""
        self.mol = Chem.MolFromSmiles("C")  # Methane
        AllChem.EmbedMolecule(self.mol, randomSeed=42)
    
    def test_identity_transformation(self):
        """Identity transformation should not change coordinates."""
        rotation = np.eye(3)
        translation = np.zeros(3)
        
        original_coords = self._get_coords(self.mol)
        transformed_mol = apply_transformation_to_molecule(self.mol, rotation, translation)
        transformed_coords = self._get_coords(transformed_mol)
        
        np.testing.assert_array_almost_equal(original_coords, transformed_coords)
    
    def test_translation_only(self):
        """Pure translation should shift all coordinates uniformly."""
        rotation = np.eye(3)
        translation = np.array([1.0, 2.0, 3.0])
        
        original_coords = self._get_coords(self.mol)
        transformed_mol = apply_transformation_to_molecule(self.mol, rotation, translation)
        transformed_coords = self._get_coords(transformed_mol)
        
        expected = original_coords + translation
        np.testing.assert_array_almost_equal(transformed_coords, expected)
    
    def test_rotation_180_z(self):
        """180° rotation around Z axis should flip X and Y."""
        rotation = np.array([
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, 1]
        ], dtype=float)
        translation = np.zeros(3)
        
        original_coords = self._get_coords(self.mol)
        transformed_mol = apply_transformation_to_molecule(self.mol, rotation, translation)
        transformed_coords = self._get_coords(transformed_mol)
        
        # X and Y should be negated, Z unchanged
        expected = original_coords.copy()
        expected[:, 0] *= -1
        expected[:, 1] *= -1
        np.testing.assert_array_almost_equal(transformed_coords, expected, decimal=5)
    
    def _get_coords(self, mol: Chem.Mol) -> np.ndarray:
        """Extract coordinates from molecule conformer."""
        conf = mol.GetConformer()
        return np.array([
            list(conf.GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ])


class TestMoleculeAlignment(unittest.TestCase):
    """Test ligand-ligand alignment functions."""
    
    def setUp(self):
        """Create test molecules."""
        self.mol_1 = Chem.MolFromSmiles("c1ccccc1")  # Benzene
        self.mol_2 = Chem.MolFromSmiles("c1ccccc1")  # Benzene (same)
        
        AllChem.EmbedMolecule(self.mol_1, randomSeed=42)
        AllChem.EmbedMolecule(self.mol_2, randomSeed=123)  # Different conformer
    
    def test_alignment_returns_scores(self):
        """Alignment should return shape similarity and RMSD."""
        shape_sim, rmsd = align_molecules(self.mol_1, self.mol_2)
        
        self.assertIsInstance(shape_sim, float)
        self.assertIsInstance(rmsd, float)
    
    def test_identical_molecule_alignment(self):
        """Aligning identical molecules should give high shape similarity."""
        mol_copy = Chem.Mol(self.mol_1)
        shape_sim, rmsd = align_molecules(self.mol_1, mol_copy)
        
        # Shape should be very high for identical molecules
        self.assertGreater(shape_sim, 0.8)


class TestFoldseekRunner(unittest.TestCase):
    """Test Foldseek integration."""
    
    def test_parse_alignment_result(self):
        """Test parsing of Foldseek output format."""
        runner = FoldseekRunner()
        
        # Create a mock output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            # Format: query,target,u,t,lddt,rmsd,alnlen
            # u = rotation matrix (9 values), t = translation (3 values)
            f.write("query.pdb\ttarget.pdb\t1,0,0,0,1,0,0,0,1\t1.0,2.0,3.0\t0.85\t1.5\t100\n")
            tmp_path = Path(f.name)
        
        try:
            result = runner._parse_alignment_result(tmp_path)
            
            self.assertIsNotNone(result)
            self.assertEqual(result.query_id, "query.pdb")
            self.assertEqual(result.target_id, "target.pdb")
            self.assertAlmostEqual(result.lddt, 0.85)
            self.assertAlmostEqual(result.rmsd, 1.5)
            self.assertEqual(result.aligned_length, 100)
            
            # Check rotation matrix (identity)
            np.testing.assert_array_almost_equal(result.rotation, np.eye(3))
            
            # Check translation
            np.testing.assert_array_almost_equal(result.translation, [1.0, 2.0, 3.0])
        finally:
            os.unlink(tmp_path)
    
    def test_parse_empty_file(self):
        """Empty output file should return None."""
        runner = FoldseekRunner()
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            tmp_path = Path(f.name)
        
        try:
            result = runner._parse_alignment_result(tmp_path)
            self.assertIsNone(result)
        finally:
            os.unlink(tmp_path)
    
    def test_parse_nonexistent_file(self):
        """Non-existent file should return None."""
        runner = FoldseekRunner()
        result = runner._parse_alignment_result(Path("/nonexistent/file.tsv"))
        self.assertIsNone(result)


class TestSystemInfo(unittest.TestCase):
    """Test SystemInfo data class and parsing."""
    
    def test_system_id_parsing(self):
        """Test that system_id is correctly parsed."""
        info = SystemInfo(
            system_id="4grb__1__1.A__1.C",
            pdb_id="4grb",
            ligand_chain="1.C",
            protein_chains=["1.A"],
        )
        
        self.assertEqual(info.system_id, "4grb__1__1.A__1.C")
        self.assertEqual(info.pdb_id, "4grb")
        self.assertEqual(info.ligand_chain, "1.C")
        self.assertEqual(info.protein_chains, ["1.A"])
    
    def test_multi_chain_system(self):
        """Test system with multiple protein chains."""
        info = SystemInfo(
            system_id="8grd__1__1.A_1.B__1.D",
            pdb_id="8grd",
            ligand_chain="1.D",
            protein_chains=["1.A", "1.B"],
        )
        
        self.assertEqual(len(info.protein_chains), 2)


class TestPLISimilarityScorer(unittest.TestCase):
    """Test the main PLISimilarityScorer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.tmp_dir = tempfile.mkdtemp()
        self.systems_dir = Path(self.tmp_dir) / "systems"
        self.systems_dir.mkdir()
        
        # Create mock annotation table
        self.annotation_table = pd.DataFrame({
            "system_id": ["test__1__1.A__1.B", "test2__1__1.A__1.B"],
            "ligand_rdkit_canonical_smiles": ["CCO", "CC(=O)O"],
        })
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
    
    def test_scorer_initialization(self):
        """Test scorer initializes correctly."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        self.assertEqual(scorer.systems_dir, self.systems_dir)
        self.assertIsNotNone(scorer.annotation_table)
    
    def test_get_system_info_from_annotation(self):
        """Test getting system info with annotation table."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        info = scorer.get_system_info("test__1__1.A__1.B")
        
        self.assertIsNotNone(info)
        self.assertEqual(info.pdb_id, "test")
        self.assertEqual(info.ligand_chain, "1.B")
        self.assertEqual(info.ligand_smiles, "CCO")
    
    def test_compute_tanimoto_from_smiles(self):
        """Test Tanimoto computation using SMILES from annotation."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        similarity = scorer.compute_tanimoto("test__1__1.A__1.B", "test2__1__1.A__1.B")
        
        # Should be computed from CCO vs CC(=O)O
        self.assertGreater(similarity, 0.0)
        self.assertLess(similarity, 1.0)
    
    def test_score_dataclass(self):
        """Test SimilarityScore dataclass."""
        score = SimilarityScore(
            query_system="sys1",
            target_system="sys2",
            tanimoto=0.75,
            sucos_protein_aligned=0.65,
            sucos_ligand_aligned=0.80,
        )
        
        self.assertEqual(score.query_system, "sys1")
        self.assertEqual(score.tanimoto, 0.75)
        self.assertTrue(np.isnan(score.foldseek_lddt))  # Default is NaN


class TestIntegrationWithMockData(unittest.TestCase):
    """Integration tests with mock PDB/SDF data."""
    
    def setUp(self):
        """Create mock system directories with test molecules."""
        self.tmp_dir = tempfile.mkdtemp()
        self.systems_dir = Path(self.tmp_dir) / "systems"
        
        # Create two mock systems
        self._create_mock_system("sys1__1__1.A__1.B", "c1ccccc1O")  # Phenol
        self._create_mock_system("sys2__1__1.A__1.B", "c1ccccc1N")  # Aniline
        
        self.annotation_table = pd.DataFrame({
            "system_id": ["sys1__1__1.A__1.B", "sys2__1__1.A__1.B"],
            "ligand_rdkit_canonical_smiles": ["c1ccccc1O", "c1ccccc1N"],
        })
    
    def tearDown(self):
        """Clean up."""
        import shutil
        shutil.rmtree(self.tmp_dir, ignore_errors=True)
    
    def _create_mock_system(self, system_id: str, smiles: str):
        """Create a mock system directory with ligand SDF."""
        system_dir = self.systems_dir / system_id
        ligand_dir = system_dir / "ligand_files"
        ligand_dir.mkdir(parents=True)
        
        # Parse system_id to get ligand chain
        parts = system_id.split("__")
        ligand_chain = parts[3]
        
        # Create molecule and save as SDF
        mol = Chem.MolFromSmiles(smiles)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        writer = Chem.SDWriter(str(ligand_dir / f"{ligand_chain}.sdf"))
        writer.write(mol)
        writer.close()
        
        # Create minimal receptor.pdb
        receptor_path = system_dir / "receptor.pdb"
        with open(receptor_path, 'w') as f:
            f.write("HEADER    TEST RECEPTOR\n")
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
            f.write("END\n")
    
    def test_load_ligand_mol(self):
        """Test loading ligand from SDF file."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        mol = scorer.load_ligand_mol("sys1__1__1.A__1.B")
        
        self.assertIsNotNone(mol)
        self.assertEqual(mol.GetNumAtoms(), 7)  # Phenol has 7 atoms
    
    def test_compute_tanimoto_from_sdf(self):
        """Test Tanimoto using loaded molecules."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        similarity = scorer.compute_tanimoto("sys1__1__1.A__1.B", "sys2__1__1.A__1.B")
        
        # Phenol and aniline have some structural similarity
        self.assertGreater(similarity, 0.0)
        self.assertLess(similarity, 1.0)
    
    def test_compute_sucos_ligand_aligned(self):
        """Test ligand-aligned SuCOS computation."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        sucos, shape, color = scorer.compute_sucos_ligand_aligned(
            "sys1__1__1.A__1.B",
            "sys2__1__1.A__1.B"
        )
        
        self.assertGreater(sucos, 0.3)  # Similar molecules
        self.assertIsInstance(shape, float)
        self.assertIsInstance(color, float)
    
    def test_score_pair_without_foldseek(self):
        """Test scoring pair without protein alignment."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        score = scorer.score_pair(
            "sys1__1__1.A__1.B",
            "sys2__1__1.A__1.B",
            compute_protein_aligned=False,  # Skip Foldseek
        )
        
        # Just check that Tanimoto was computed (not NaN)
        self.assertGreater(score.tanimoto, 0.0)
        self.assertGreater(score.sucos_ligand_aligned, 0.3)
        self.assertTrue(np.isnan(score.sucos_protein_aligned))  # Not computed
    
    def test_score_batch(self):
        """Test batch scoring."""
        scorer = PLISimilarityScorer(
            systems_dir=self.systems_dir,
            annotation_table=self.annotation_table,
        )
        
        pairs = [
            ("sys1__1__1.A__1.B", "sys2__1__1.A__1.B"),
            ("sys2__1__1.A__1.B", "sys1__1__1.A__1.B"),  # Reverse
        ]
        
        df = scorer.score_batch(pairs, compute_protein_aligned=False, progress=False)
        
        self.assertEqual(len(df), 2)
        self.assertIn("tanimoto", df.columns)
        self.assertIn("sucos_ligand_aligned", df.columns)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""
    
    def test_nonexistent_system(self):
        """Non-existent system should return None/NaN gracefully."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            scorer = PLISimilarityScorer(
                systems_dir=Path(tmp_dir),
                annotation_table=pd.DataFrame({"system_id": [], "ligand_rdkit_canonical_smiles": []}),
            )
            
            mol = scorer.load_ligand_mol("nonexistent__1__1.A__1.B")
            self.assertIsNone(mol)
    
    def test_empty_annotation_table(self):
        """Empty annotation table should be handled."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            scorer = PLISimilarityScorer(
                systems_dir=Path(tmp_dir),
                annotation_table=pd.DataFrame(),
            )
            
            info = scorer.get_system_info("test__1__1.A__1.B")
            self.assertIsNotNone(info)
            self.assertIsNone(info.ligand_smiles)
    
    def test_invalid_smiles_in_annotation(self):
        """Invalid SMILES in annotation should be handled."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            annotation = pd.DataFrame({
                "system_id": ["test__1__1.A__1.B"],
                "ligand_rdkit_canonical_smiles": ["invalid_smiles"],
            })
            
            scorer = PLISimilarityScorer(
                systems_dir=Path(tmp_dir),
                annotation_table=annotation,
            )
            
            fp = scorer.get_fingerprint("test__1__1.A__1.B")
            self.assertIsNone(fp)


if __name__ == "__main__":
    unittest.main(verbosity=2)
