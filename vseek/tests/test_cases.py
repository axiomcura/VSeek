import unittest
import numpy as np
import pandas as pd
import json
from pathlib import Path

# VSeek imports
import vseek.common.io_files as vfiles
import vseek.common.loader as vloader
import vseek.common.vseek_paths as vsp
from vseek.utils.sequence_io import SequenceIO
from vseek.utils.sra_callers import download_fasta
from vseek.utils.vseek_analysis import dynamic_hamming
from vseek.utils.vseek_plots import plot_viral_composition, bat_country_geo_plot
from vseek.apis.ncbi import get_all_viral_accessions, get_viral_genes, get_viral_genomes
from scipy.spatial.distance import hamming
from vseek.utils.vseek_analysis import dynamic_hamming, hamming_distance_score


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

    def test_hamming_score(self):
        seq1 = "AAATTCCC"
        seq2 = "ACATGCCC"
        expected_resp = 2 / 8
        test_resp = hamming_distance_score(seq1, seq2)

        self.assertEqual(expected_resp, test_resp)

    def test_loading_viral_genes(self):
        exp_gene = "TGGCTACTTGGACCCCTAACACTGGACGGCTTTATCTTCCTCCTGCTAAACCTGTGGCGACTGTTCTATCGACTGATGATTATATTGTACCGACGAATCT"
        test_gene = vloader.load_viral_genes("NC_013035")
        self.assertEqual(exp_gene, test_gene[0][:100])


class TestAnalysis(unittest.TestCase):
    def test_hamming_calculation(self):
        seq1 = "TTAAA"
        seq2 = "TTACG"
        seq1_arr = [c for c in seq1]
        seq2_arr = [c for c in seq2]

        scipy_hamming = hamming(seq1_arr, seq2_arr)
        test_hamming = hamming_distance_score(seq1, seq2)

        self.assertEqual(scipy_hamming, test_hamming)

    def test_dynamic_hamming(self):

        ref_seq = "TGCTATTCAGGGCTTGACCAACACTGGATTGCTTTTCACTTAAAGTATTATGCACGACAGGGTGCGTGTACCATGTAAACCTGTTATAACTTACCTCAGA"
        read = "TTAAA"

        # using scipy as control
        start_idx = 0
        end_idx = len(read)
        scores = []
        while end_idx - 1 < len(ref_seq):
            ref_piece_arr = [c for c in ref_seq[start_idx:end_idx]]
            read_arr = [c for c in read]
            score = hamming(read_arr, ref_piece_arr)

            # increment
            start_idx += 1
            end_idx += 1
            scores.append(score)

        expected_final_score = max(scores)
        test_final_score = dynamic_hamming(read, ref_seq)

        self.assertEqual(1.0 - expected_final_score, test_final_score)


class TestProfilerAndLoader(unittest.TestCase):
    def test_viral_accesion(self):
        def rel_abundance(counts, count_sum):
            return round((counts / count_sum))

        rel_threshold = 1.0
        # loading viral counts
        ncbi_virus_data = vloader.load_bat_virus_data()
        with open("./SRR12464727_viral_composition_counts.json") as outfile:
            viral_count_data = json.load(outfile)
        viral_count_df = pd.DataFrame.from_dict(
            data=viral_count_data, orient="index"
        ).reset_index()
        viral_count_df.columns = ["accession", "counts"]
        viral_count_sum = np.sum(viral_count_df["counts"])
        viral_count_df["rel_abundance"] = viral_count_df["counts"].apply(
            lambda count: rel_abundance(counts=count, count_sum=viral_count_sum)
        )
        viral_count_df = viral_count_df.sort_values(by="rel_abundance", ascending=False)

        selected_counts = viral_count_df.loc[
            viral_count_df["rel_abundance"] >= rel_threshold
        ]
        viral_taxa = ncbi_virus_data[["Representative", "family", "genus"]]
        viral_taxa.columns = ["accession", "family", "genus"]
        identified_virus_df = viral_taxa.merge(selected_counts, on="accession")
        identified_virus_df = identified_virus_df.sort_values(
            by="rel_abundance", ascending=False
        )
        identified_virus_df = identified_virus_df.set_index("accession")
        exp_pout = ["Herpesviridae", "Cytomegalovirus", 1591, 1]
        test_out = identified_virus_df.iloc[0].values.tolist()

        self.assertEqual(exp_pout, test_out)

    def test_ppi_profile(self):
        test_ds = [
            [
                9606,
                "ENSP00000000233",
                10335,
                "NP04_VZVD",
                437,
                "Varicella-zoster virus (strain Dumas)",
            ],
            [
                9606,
                "ENSP00000000233",
                10255,
                "F11_VAR67",
                188,
                "Variola virus (isolate Human/India/Ind3/1967)",
            ],
            [
                9606,
                "ENSP00000000233",
                10335,
                "GM_VZVD",
                336,
                "Varicella-zoster virus (strain Dumas)",
            ],
        ]
        # loader contains internal functions:
        ppi_df = vloader.load_human_ppi()

        self.assertEqual(test_ds, ppi_df.iloc[:3].values.tolist())

    def test_Dbatvir_vir_(self):
        test_data = [
            [
                "Bat adeno-associated virus 07YN",
                "Parvoviridae",
                "unclassified Chiroptera",
                np.nan,
                np.nan,
                "Feces",
                "2007",
                "China",
                "Cap",
                "Unpublished",
            ],
            [
                "Bat adeno-associated virus 09YN",
                "Parvoviridae",
                "unclassified Chiroptera",
                np.nan,
                np.nan,
                "Feces",
                "2009",
                "China",
                "Cap",
                "Unpublished",
            ],
        ]

        exp = vloader.load_dbat_vir_db()[:2].values.tolist()
        self.assertEqual(exp, test_data)

    def test_location_cods(self):

        test_set = [
            ["Afghanistan", "AF", "AFG"],
            ["Åland Islands", "AX", "ALA"],
            ["Albania", "AL", "ALB"],
            ["Algeria", "DZ", "DZA"],
            ["American Samoa", "AS", "ASM"],
        ]

        exp = vloader.load_iban_iso_codes().iloc[:5].values.tolist()
        self.assertEqual(test_set, exp)

    def testing_bat_country_profile(self):
        test_viral_loc = [
            [
                "Spain",
                "Herpesviridae",
                "Eptesicus isabellinus Hypsugo savii Myotis alcathoe Myotis bechsteinii Myotis blythii Myotis capaccinii Myotis daubentonii Myotis emarginatus Myotis escalerai Myotis myotis Myotis mystacinus Nyctalus lasiopterus Nyctalus leisleri Pipistrellus kuhlii Pipistrellus pygmaeus Plecotus austriacus Rhinolophus hipposideros Tadarida teniotis",
                18,
            ]
        ]

        def rel_abundance(counts, counts_sum):
            return round((counts / counts_sum))

        rel_threshold = 1.0
        # loading viral counts
        ncbi_virus_data = vloader.load_bat_virus_data()
        with open("./SRR12464727_viral_composition_counts.json") as outfile:
            viral_count_data = json.load(outfile)

        viral_count_df = pd.DataFrame.from_dict(
            data=viral_count_data, orient="index"
        ).reset_index()
        viral_count_df.columns = ["accession", "counts"]
        viral_count_sum = np.sum(viral_count_df["counts"])
        viral_count_df["rel_abundance"] = viral_count_df["counts"].apply(
            lambda count: rel_abundance(counts=count, counts_sum=viral_count_sum)
        )
        viral_count_df = viral_count_df.sort_values(by="rel_abundance", ascending=False)

        selected_counts = viral_count_df.loc[
            viral_count_df["rel_abundance"] >= rel_threshold
        ]
        viral_taxa = ncbi_virus_data[["Representative", "family", "genus"]]
        viral_taxa.columns = ["accession", "family", "genus"]
        identified_virus_df = viral_taxa.merge(selected_counts, on="accession")
        identified_virus_df = identified_virus_df.sort_values(
            by="rel_abundance", ascending=False
        )
        identified_virus_df = identified_virus_df.set_index("accession")

        geolocations_df = vloader.load_geolocations()
        dbat_vir_df = vloader.load_dbat_vir_db()
        found_viral_fam = identified_virus_df["family"].unique().tolist()
        associated_bats_df = dbat_vir_df.loc[
            dbat_vir_df["Viral family"].isin(found_viral_fam)
        ]

        # removing any duplicates
        group = associated_bats_df.groupby("From Bat")

        list_dfs = []
        for name, df in group:
            new_df = df.drop_duplicates(subset=["Viral family", "From Bat"])
            list_dfs.append(new_df)

        associated_bats_df = pd.concat(list_dfs)
        associated_bats_df["Sampling country"].unique()

        # grouping by country and condensing data
        country_group = associated_bats_df.groupby("Sampling country")
        bat_country_data = []
        for name, df in country_group:

            if "congo" in name.lower():
                name = "Congo [Republic]"
            elif "united states of america" in name.lower():
                name = "United States"

            all_bats_str = " ".join(df["From Bat"].tolist())
            counts = len(df["From Bat"].tolist())
            viral_fam = ", ".join(df["Viral family"].unique().tolist())
            result = [name, viral_fam, all_bats_str, counts]

            bat_country_data.append(result)

        # constructing geo dataframe
        all_bats_country_df = pd.DataFrame(
            data=bat_country_data,
            columns=["country", "viral_fam", "all_bats", "bat_counts"],
        )
        all_bats_country_df = all_bats_country_df.sort_values(
            by="bat_counts", ascending=False
        )

        test_loc_profile = all_bats_country_df.iloc[:1].values.tolist()
        self.assertEqual(test_viral_loc, test_loc_profile)
