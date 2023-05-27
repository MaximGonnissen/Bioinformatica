import unittest

from blosum import BLOSUM

from src.fasta_parser.fasta_parser import parse
from src.psa.smith_waterman import SmithWatermanPSASolver


class TestSmithWaterman(unittest.TestCase):
    """
    Tests for the Smith-Waterman algorithm.
    """

    def test_smith_waterman(self):
        """
        Test the Smith-Waterman algorithm.
        """
        sequences = parse("psa_tests/test_inputs/test.fasta")

        sequence_ids = list(sequences.keys())
        sequence2 = sequences[sequence_ids[1]]
        sequence1 = sequences[sequence_ids[0]]

        config = {
            'indel': 2
        }

        solver = SmithWatermanPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, 7)
        self.assertEqual(len(alignments), 1)
        self.assertEqual(alignments[0][1], "TEAFF")
        self.assertEqual(alignments[0][0], "SKIIF")

    def test_smith_waterman_2(self):
        """
        Test the Smith-Waterman algorithm with an indel penalty of 1.
        """
        sequences = parse("psa_tests/test_inputs/test.fasta")

        sequence_ids = list(sequences.keys())
        sequence2 = sequences[sequence_ids[1]]
        sequence1 = sequences[sequence_ids[0]]

        config = {
            'indel': 1
        }

        solver = SmithWatermanPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, 8)
        self.assertEqual(len(alignments), 4)

        correct_pairs = [
            ("SSASKIIF", "TEA--F-F"),
            ("SSASKIIF", "TEA---FF"),
            ("SS-ASKIIF", "NTEA---FF"),
            ("SS-ASKIIF", "NTEA--F-F")
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)
