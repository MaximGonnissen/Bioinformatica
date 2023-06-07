import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.scoring_matrix.scoring_matrix import ScoringMatrix
from src.msa.smith_waterman import SmithWatermanMSASolver


class TestSmithWatermanMSA(unittest.TestCase):
    """
    Tests for the Smith-Waterman multiple sequence alignment solver.
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
        self.scoring_matrix = ScoringMatrix(self.sequence_values)

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
        solver = SmithWatermanMSASolver(self.config)
        score, alignments = solver.solve(self.sequence_values)

        # Correct alignments and score unknown.
        # TODO: Find correct score.


class TestSmithWatermanMSA2D(TestSmithWatermanMSA):
    """
    Tests for the Smith-Waterman multiple sequence alignment solver with 2D input.
    """
    test_input = "msa_tests/test_inputs/test_2d.fasta"

    def test_solve(self) -> None:
        """
        Test the solve method.
        """
        solver = SmithWatermanMSASolver(self.config)
        score, alignments = solver.solve(self.sequence_values)

        correct_alignments = [
            ("FGSGTRL", "FGQGTRL")
        ]

        self.assertEqual(score, 28)
        self.check_alignments(alignments, correct_alignments)
