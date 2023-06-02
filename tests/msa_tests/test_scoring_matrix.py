import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.scoring_matrix.scoring_matrix import ScoringMatrix


class TestScoringMatrix(unittest.TestCase):
    """
    Tests for the MSA scoring matrix.
    """

    def test_initialise(self):
        """
        Test initialisation of the scoring matrix.
        """
        sequences = parse("msa_tests/test_inputs/test.fasta")

        sequences = [sequence for sequence in sequences.values()]

        scoring_matrix = ScoringMatrix(sequences)

        self.assertEqual(scoring_matrix.shape, (21, 16, 16))

        for i in range(21):
            for j in range(16):
                for k in range(16):
                    self.assertEqual(scoring_matrix[i, j, k], 0)
