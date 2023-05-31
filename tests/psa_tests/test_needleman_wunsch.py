import unittest

from blosum import BLOSUM

from src.fasta_parser.fasta_parser import parse
from src.psa.needleman_wunsch import NeedlemanWunschPSASolver


class TestNeedlemanWunsch(unittest.TestCase):
    """
    Tests for the Needleman-Wunsch algorithm.
    """

    def test_needleman_wunsch(self):
        """
        Test the Needleman-Wunsch algorithm.
        """
        sequences = parse("psa_tests/test_inputs/test.fasta")

        sequence_ids = list(sequences.keys())
        sequence2 = sequences[sequence_ids[1]]
        sequence1 = sequences[sequence_ids[0]]

        config = {
            'indel': 2
        }

        solver = NeedlemanWunschPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, -5)
        self.assertEqual(len(alignments), 2)

        correct_pairs = [
            ("GYSSASKIIF----", "N-TEA---FFGQGT"),
            ("GYSSASKIIF----", "N-TEA--F-FGQGT"),
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)
