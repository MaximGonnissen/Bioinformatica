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
            'indel': -2
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
            'indel': -1
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

    def test_psa_init_match_mismatch(self):
        """
        Test that a PSA solver's score matrix uses match and mismatch values if they are provided.
        """
        config = {
            'match': 5,
            'mismatch': -5,
            'indel': -2
        }

        solver = SmithWatermanPSASolver(config=config)

        self.assertEqual(solver.substitution_matrix['A']['A'], 5)
        self.assertEqual(solver.substitution_matrix['A']['C'], -5)

        for key in solver.substitution_matrix.keys():
            for other_key in solver.substitution_matrix[key].keys():
                if key != other_key:
                    self.assertEqual(solver.substitution_matrix[key][other_key], -5)
                else:
                    self.assertEqual(solver.substitution_matrix[key][other_key], 5)

    def test_smith_waterman_assignment(self):
        """
        Test the Smith-Waterman algorithm with the assignment example.
        """
        sequences = parse("psa_tests/test_inputs/assignment_test.fasta")

        sequences = list(sequences.values())

        sequence1 = sequences[0]
        sequence2 = sequences[1]

        config = {
            "match": 5,
            "mismatch": -2,
            "indel": -4,
            "two gaps": 0
        }

        solver = SmithWatermanPSASolver(config=config)

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, 28)
        self.assertEqual(len(alignments), 1)

        correct_pairs = [
            ("FGSGTRL", "FGQGTRL")
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)
