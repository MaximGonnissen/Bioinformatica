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
            "indel": -4,
        }

        solver = NeedlemanWunschPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, -15)
        self.assertEqual(len(alignments), 3)

        correct_pairs = [
            ("GYSS-ASKIIF", "NTEAFFGQGT-"),
            ("GYSSA-SKIIF", "NTEAFFGQGT-"),
            ("GYSSA--SKIIF", "N-TEAFFGQGT-")
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)

    def test_needleman_wunsch_2(self):
        """
        Test the Needleman-Wunsch algorithm with an indel penalty of 3.
        """
        sequences = parse("psa_tests/test_inputs/test.fasta")

        sequence_ids = list(sequences.keys())
        sequence2 = sequences[sequence_ids[1]]
        sequence1 = sequences[sequence_ids[0]]

        config = {
            'indel': -3
        }

        solver = NeedlemanWunschPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, -11)
        self.assertEqual(len(alignments), 1)

        correct_pairs = [
            ("GYSSA--SKIIF", "N-TEAFFGQGT-")
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)

    def test_needleman_wunsch_3(self):
        """
        Test the Needleman-Wunsch algorithm with an indel penalty of 15.
        """
        sequences = parse("psa_tests/test_inputs/test.fasta")

        sequence_ids = list(sequences.keys())
        sequence2 = sequences[sequence_ids[1]]
        sequence1 = sequences[sequence_ids[0]]

        config = {
            'indel': -15
        }

        solver = NeedlemanWunschPSASolver(config=config, substitution_matrix=BLOSUM(62))

        score, alignments = solver.solve(sequence_1=sequence1, sequence_2=sequence2)

        self.assertEqual(score, -16)
        self.assertEqual(len(alignments), 1)

        correct_pairs = [
            ("GYSSASKIIF", "NTEAFFGQGT")
        ]

        for alignment in alignments:
            self.assertIn(alignment, correct_pairs)
