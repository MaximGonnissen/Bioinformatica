import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.scoring_matrix.scoring_matrix import ScoringMatrix


class TestSmithWatermanMSA(unittest.TestCase):
    """
    Tests for the Smith-Waterman multiple sequence alignment solver.
    """

    def setUp(self) -> None:
        """
        Set up the test case.
        """
        self.sequences = parse("msa_tests/test_inputs/test.fasta")
        self.sequence_values = list(self.sequences.values())
        self.scoring_matrix = ScoringMatrix(self.sequence_values)
        self.config = {
            "indel": 4,
            "two gaps": 0
        }

    def test_solve(self):
        """
        Test the Smith-Waterman multiple sequence alignment solver.
        """
        pass
