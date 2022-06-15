import unittest

from scipy.spatial.distance import hamming
from vseek.utils.vseek_analysis import dynamic_hamming, hamming_distance_score

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

