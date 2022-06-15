import unittest
from pathlib import Path

#VSeek Imports
import vseek.common.vseek_paths as vsp

class TestFunction(unittest.TestCase):
    def test_directories(self):

        # expected
        root_path = "./VSeek"
        exp_root_path = str(Path("./").absolute())
        exp_db_dir = str(Path("./db").absolute())
        exp_fasta_dir = str(Path("./results/fasta_files").absolute())
        exp_genome_dir = str(Path("./db/genome").absolute())
        exp_ppi_dir = str(Path("./db/protein_interactions").absolute())
        exp_results_dir = str(Path("./results").absolute())
        exp_sra_dir = str(Path("./results/SRA_prefetch").absolute())

        # testing paths
        self.assertEqual(root_path, vsp.relative_root_path())
        self.assertEqual(exp_db_dir, vsp.db_path())
        self.assertEqual(exp_genome_dir, vsp.genome_db_path())
        self.assertEqual(exp_results_dir, vsp.results_dir())
        self.assertEqual(exp_fasta_dir, vsp.metagenome_path())
        self.assertEqual(exp_ppi_dir, vsp.ppi_db_path())
        self.assertEqual(exp_sra_dir, vsp.prefetch_path())

        



