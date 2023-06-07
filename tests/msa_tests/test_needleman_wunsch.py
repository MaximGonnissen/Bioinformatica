import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.needleman_wunsch import NeedlemanWunschMSASolver
from src.msa.scoring_matrix.scoring_matrix import ScoringMatrix


class TestNeedlemanWunschMSA(unittest.TestCase):
    """
    Tests for the Needleman-Wunsch multiple sequence alignment solver.
    """
    test_input = "msa_tests/test_inputs/test.fasta"
    config = {
        "match": 5,
        "mismatch": -2,
        "indel": -4,
        "two gaps": 0
    }

    def setUp(self) -> None:
        """
        Set up the test case.
        """
        self.sequences = parse(self.test_input)
        self.sequence_values = list(self.sequences.values())

    def check_alignments(self, alignments: list[tuple[str, ...]], correct_alignments: list[tuple[str, ...]]) -> None:
        """
        Check that the alignments are correct.
        :param alignments: The alignments to check.
        :param correct_alignments: The correct alignments.
        """
        for alignment in correct_alignments:
            self.assertIn(alignment, alignments)

    def test_solve(self) -> None:
        """
        Test the solve method.
        """
        solver = NeedlemanWunschMSASolver(self.config)
        score, alignments = solver.solve(self.sequence_values)


class TestNeedlemanWunschMSA2D(TestNeedlemanWunschMSA):
    """
    Tests for the Smith-Waterman multiple sequence alignment solver with 2D input.
    """
    test_input = "msa_tests/test_inputs/test_2d.fasta"

    def test_solve(self) -> None:
        """
        Test the solve method.

        Note: This alignment seems to have over 50 possible alignments, so this was purely used for score testing.
        """
        solver = NeedlemanWunschMSASolver(self.config)
        score, alignments = solver.solve(self.sequence_values)

        self.assertEqual(score, -1)

    def test_solve_short(self):
        """
        Test the solve method with a short input.
        """
        sequences = [
            "GYSSA",
            "NTEAFF"
        ]

        solver = NeedlemanWunschMSASolver(self.config)
        score, alignments = solver.solve(sequences)

        correct_alignments = [
            ("GYSSA--", "-NTEAFF"),
            ("GYSSA--", "N-TEAFF"),
            ("GYSSA--", "NT-EAFF"),
            ("GYSSA--", "NTE-AFF")
        ]

        self.assertEqual(score, -13)

        self.check_alignments(alignments, correct_alignments)

    def test_solve_short2(self):
        """
        Test the solve method with another short input, with only a single possible alignment.
        """
        sequences = [
            "AATCG",
            "AACG",
        ]

        solver = NeedlemanWunschMSASolver(self.config)
        score, alignments = solver.solve(sequences)

        correct_alignments = [
            ("AATCG", "AA-CG")
        ]

        self.assertEqual(score, 16)

        self.check_alignments(alignments, correct_alignments)

    def test_solve_short3(self):
        """
        Test the solve method with another short input, with two possible alignments.
        """
        sequences = [
            "AATCGC",
            "AACGAA",
        ]

        solver = NeedlemanWunschMSASolver(self.config)
        score, alignments = solver.solve(sequences)

        correct_alignments = [
            ("AATCG-C", "AA-CGAA"),
            ("AATCGC-", "AA-CGAA")
        ]

        self.assertEqual(score, 10)

        self.check_alignments(alignments, correct_alignments)
