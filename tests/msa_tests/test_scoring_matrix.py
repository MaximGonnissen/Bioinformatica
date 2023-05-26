import unittest

from src.msa.enums import Direction
from src.msa.scoring_matrix import ScoringMatrix


class TestScoringMatrix(unittest.TestCase):
    """
    Tests for the ScoringMatrix class.
    """

    def setUp(self) -> None:
        """
        Set up the test case.
        """
        self.matrix = ScoringMatrix("ABC", "DEF")

    def test_get_score(self):
        """
        Test the get_score method.
        """
        for i in range(self.matrix.width):
            for j in range(self.matrix.height):
                self.assertEqual(self.matrix.get_score(i, j), 0)

    def test_get_traceback(self):
        """
        Test the get_traceback method.
        """
        for i in range(self.matrix.width):
            for j in range(self.matrix.height):
                self.assertEqual(self.matrix.get_traceback(i, j), [])

    def test_set_score(self):
        """
        Test the set_score method.
        """
        for i in range(self.matrix.width):
            for j in range(self.matrix.height):
                self.matrix.set_score(i, j, 1)
                self.assertEqual(self.matrix.get_score(i, j), 1)

    def test_set_traceback(self):
        """
        Test the set_traceback method.
        """
        for i in range(self.matrix.width):
            for j in range(self.matrix.height):
                self.matrix.set_traceback(i, j, [Direction.DIAGONAL])
                self.assertEqual(self.matrix.get_traceback(i, j), [Direction.DIAGONAL])

    def test_initialise(self):
        """
        Test the initialise method.
        """
        self.matrix.initialise()
        for i in range(self.matrix.width):
            for j in range(self.matrix.height):
                self.assertEqual(self.matrix.get_score(i, j), 0)
                self.assertEqual(self.matrix.get_traceback(i, j), [])

    def test_get_max_score(self):
        """
        Test the get_max_score method.
        """
        self.assertEqual(self.matrix.max_score(), 0)
        self.assertEqual(self.matrix.max_score_index_multiple(),
                         [(x, y) for x in range(self.matrix.width) for y in range(self.matrix.height)])

        self.matrix.set_score(0, 0, 1)
        self.assertEqual(self.matrix.max_score(), 1)
        self.assertEqual(self.matrix.max_score_index(), (0, 0))

        self.matrix.set_score(0, 1, 2)
        self.assertEqual(self.matrix.max_score(), 2)
        self.assertEqual(self.matrix.max_score_index(), (0, 1))

        self.matrix.set_score(1, 0, 3)
        self.assertEqual(self.matrix.max_score(), 3)
        self.assertEqual(self.matrix.max_score_index(), (1, 0))
