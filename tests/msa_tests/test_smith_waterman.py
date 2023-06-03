import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.smith_waterman import SmithWatermanMSASolver


class TestSmithWatermanMSA(unittest.TestCase):
    """
    Tests for the Smith-Waterman multiple sequence alignment solver.
    """

    def test_solve(self):
        """
        Test the Smith-Waterman multiple sequence alignment solver.
        """
        sequences = parse("msa_tests/test_inputs/test.fasta")

        sequences = [sequence for sequence in sequences.values()]

        config = {
            "indel": 4,
            "two gaps": 0
        }

        solver = SmithWatermanMSASolver(config)

        alignments = solver.solve(sequences)

        # TODO: Write actual asserts -- need correct data to use as verification data.

        # print(alignments)
