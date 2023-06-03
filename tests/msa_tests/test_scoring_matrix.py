import unittest

from src.fasta_parser.fasta_parser import parse
from src.msa.scoring_matrix.scoring_matrix import ScoringMatrix, ScoringMatrixEntry


class TestScoringMatrix(unittest.TestCase):
    """
    Tests for the MSA scoring matrix.
    """

    def setUp(self) -> None:
        """
        Set up the test case.
        """
        self.sequences = parse("msa_tests/test_inputs/test.fasta")
        self.sequence_values = [sequence for sequence in self.sequences.values()]
        self.scoring_matrix = ScoringMatrix(self.sequence_values)

    def test_initialise(self):
        """
        Test initialisation of the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.shape, (21, 16, 16))

        for i in range(21):
            for j in range(16):
                for k in range(16):
                    self.assertEqual(self.scoring_matrix[i, j, k].score, 0)

    def test_set_direct(self):
        """
        Test the setter for the scoring matrix.
        """
        self.scoring_matrix[0, 0, 0] = 1
        self.assertEqual(self.scoring_matrix[0, 0, 0].score, 1)

        self.scoring_matrix[0, 0, 1] = 2
        self.assertEqual(self.scoring_matrix[0, 0, 0].score, 1)
        self.assertEqual(self.scoring_matrix[0, 0, 1].score, 2)

        self.scoring_matrix[0, 0, 2] = 3
        self.assertEqual(self.scoring_matrix[0, 0, 0].score, 1)
        self.assertEqual(self.scoring_matrix[0, 0, 1].score, 2)
        self.assertEqual(self.scoring_matrix[0, 0, 2].score, 3)

        self.scoring_matrix[0, 0, 0] = ScoringMatrixEntry(4, [1, 2, 3])
        self.assertEqual(self.scoring_matrix[0, 0, 0].score, 4)
        self.assertEqual(self.scoring_matrix[0, 0, 0].traceback, [1, 2, 3])
        self.assertEqual(self.scoring_matrix[0, 0, 1].score, 2)
        self.assertEqual(self.scoring_matrix[0, 0, 2].score, 3)

        self.assertNotEqual(self.scoring_matrix[0, 0, 0].score, 1)

    def test_set_score(self):
        """
        Test setting a score in the scoring matrix.
        """
        self.scoring_matrix.set_score(0, 0, 0, score=1)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)

        self.scoring_matrix.set_score(0, 0, 1, score=2)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 2)

        self.scoring_matrix.set_score(0, 0, 2, score=3)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 2)
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 2), 3)

    def test_add_traceback(self):
        """
        Test adding a traceback to the scoring matrix.
        """
        self.scoring_matrix.add_traceback(0, 0, 0, traceback_direction=(1, 1, 1))
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.scoring_matrix.add_traceback(0, 0, 0, traceback_direction=(2, 2, 2))
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1), (2, 2, 2)])

    def test_max_score(self):
        """
        Test finding the maximum score in the scoring matrix.
        """
        self.scoring_matrix.set_score(0, 0, 0, score=1)
        self.scoring_matrix.set_score(0, 0, 1, score=2)
        self.scoring_matrix.set_score(0, 1, 0, score=3)
        self.scoring_matrix.set_score(0, 1, 1, score=4)

        self.assertEqual(self.scoring_matrix.get_max_score(), 4)

    def test_max_score_index(self):
        """
        Test finding the index of the maximum score in the scoring matrix.
        """
        self.scoring_matrix.set_score(0, 0, 0, score=1)
        self.scoring_matrix.set_score(0, 0, 1, score=2)
        self.scoring_matrix.set_score(0, 1, 0, score=3)
        self.scoring_matrix.set_score(0, 1, 1, score=4)

        self.assertEqual(self.scoring_matrix.get_max_index(), (0, 1, 1))

        self.scoring_matrix.set_score(3, 1, 2, score=5)

        self.assertEqual(self.scoring_matrix.get_max_index(), (3, 1, 2))

    def test_get_traceback(self):
        """
        Test getting the traceback from the scoring matrix.
        """
        self.scoring_matrix.add_traceback(0, 0, 0, traceback_direction=(1, 1, 1))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.scoring_matrix.add_traceback(0, 0, 0, traceback_direction=(2, 2, 2))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1), (2, 2, 2)])

    def test_mixed_set(self):
        """
        Test that setting both scores and tracebacks don't cause issues.
        """
        self.scoring_matrix.set_score(0, 0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [])

        self.scoring_matrix.add_traceback(0, 0, 0, traceback_direction=(1, 1, 1))

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 1), [])

        self.scoring_matrix.set_score(0, 0, 1, score=2)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 1), [])

        self.scoring_matrix.add_traceback(0, 0, 1, traceback_direction=(2, 2, 2))

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 1), [(2, 2, 2)])

        self.scoring_matrix.set_traceback(0, 0, 1, traceback=[(3, 3, 3)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 1), [(3, 3, 3)])

    def test_invalid_direct_set(self):
        """
        Test that setting invalid values directly to the matrix raises an error.
        """
        with self.assertRaises(TypeError):
            self.scoring_matrix[0, 0, 0] = "A"

        with self.assertRaises(TypeError):
            self.scoring_matrix[0, 0, 0] = object

        with self.assertRaises(TypeError):
            self.scoring_matrix[0, 0, 0] = None

    def test_direct_set_traceback(self):
        """
        Test that setting tracebacks directly to the matrix works.
        """
        self.scoring_matrix[0, 0, 0] = [(1, 1, 1)]

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(1, 1, 1)])

        self.scoring_matrix[0, 0, 0] = [(2, 2, 2)]

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(2, 2, 2)])

        self.scoring_matrix[0, 0, 0] = [(3, 3, 3)]

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [(3, 3, 3)])

    def test_get_corner_index(self):
        """
        Test that getting the corner index works.
        """
        matrix_shape = self.scoring_matrix.shape

        self.assertEqual(self.scoring_matrix.get_corner_index(),
                         (matrix_shape[0] - 1, matrix_shape[1] - 1, matrix_shape[2] - 1))

    def test_needleman_wunsch_init(self):
        """
        Test that the Needleman-Wunsch initialization works in 3D.
        """
        gap_penalty = 1

        self.scoring_matrix.init_needleman_wunsch(gap_penalty=gap_penalty)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0), [])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 1), -1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 1), [(0, 0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(0, 1, 0), -1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 1, 0), [(0, 0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(1, 0, 0), -1)
        self.assertEqual(self.scoring_matrix.get_traceback(1, 0, 0), [(0, 0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(1, 1, 1), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(1, 1, 1), [])

        self.assertEqual(self.scoring_matrix.get_score(2, 0, 0), -2)
        self.assertEqual(self.scoring_matrix.get_traceback(2, 0, 0), [(1, 0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(0, 2, 0), -2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 2, 0), [(0, 1, 0)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 2), -2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 2), [(0, 0, 1)])

        self.assertEqual(self.scoring_matrix.get_score(2, 2, 2), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(2, 2, 2), [])

        self.assertEqual(self.scoring_matrix.get_score(3, 0, 0), -3)
        self.assertEqual(self.scoring_matrix.get_traceback(3, 0, 0), [(2, 0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(0, 3, 0), -3)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 3, 0), [(0, 2, 0)])

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 3), -3)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 3), [(0, 0, 2)])


class TestScoringMatrix2D(unittest.TestCase):
    """
    Test an MSA scoring matrix for 2 sequences.
    """

    def setUp(self):
        self.sequences = ["ABC", "ABD"]
        self.scoring_matrix = ScoringMatrix(self.sequences)

    def test_get_score(self):
        """
        Test getting a score from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_score(0, 0), 0)

        self.scoring_matrix.set_score(0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 1)

    def test_get_max_score(self):
        """
        Test getting the maximum score from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_max_score(), 0)

        self.scoring_matrix.set_score(0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_max_score(), 1)

        self.scoring_matrix.set_score(1, 1, score=2)

        self.assertEqual(self.scoring_matrix.get_max_score(), 2)

    def test_get_max_index(self):
        """
        Test getting the maximum index from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_max_index(), (0, 0))

        self.scoring_matrix.set_score(1, 1, score=2)

        self.assertEqual(self.scoring_matrix.get_max_index(), (1, 1))

    def test_get_traceback(self):
        """
        Test getting the traceback from the scoring matrix.
        """
        self.scoring_matrix.add_traceback(0, 0, traceback_direction=(1, 1))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [(1, 1)])

        self.scoring_matrix.add_traceback(0, 0, traceback_direction=(2, 2))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [(1, 1), (2, 2)])

    def test_mixed_set(self):
        """
        Test that setting both scores and tracebacks don't cause issues.
        """
        self.scoring_matrix.set_score(0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [])

        self.scoring_matrix.add_traceback(0, 0, traceback_direction=(1, 1))

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [(1, 1)])

        self.scoring_matrix.set_score(0, 0, score=2)

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [(1, 1)])

    def test_direct_set(self):
        """
        Test that setting scores directly to the matrix works.
        """
        self.scoring_matrix[0, 0] = 1

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [])
        self.assertEqual(self.scoring_matrix[0, 0], 1)
        self.assertEqual(self.scoring_matrix[0, 1], 0)

        self.scoring_matrix[0, 0] = 2

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [])
        self.assertEqual(self.scoring_matrix[0, 0], 2)
        self.assertEqual(self.scoring_matrix[0, 1], 0)

        self.scoring_matrix[0, 0] = [1, 2, 3]

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [1, 2, 3])
        self.assertEqual(self.scoring_matrix[0, 0], 2)
        self.assertEqual(self.scoring_matrix[0, 1], 0)

    def test_needleman_wunsch_init(self):
        """
        Test that the Needleman-Wunsch initialization works in 2D.
        """
        gap_penalty = 1

        self.scoring_matrix.init_needleman_wunsch(gap_penalty=gap_penalty)

        self.assertEqual(self.scoring_matrix.get_score(0, 0), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0), [])

        self.assertEqual(self.scoring_matrix.get_score(0, 1), -1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 1), [(0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(1, 0), -1)
        self.assertEqual(self.scoring_matrix.get_traceback(1, 0), [(0, 0)])

        self.assertEqual(self.scoring_matrix.get_score(1, 1), 0)
        self.assertEqual(self.scoring_matrix.get_traceback(1, 1), [])

        self.assertEqual(self.scoring_matrix.get_score(2, 0), -2)
        self.assertEqual(self.scoring_matrix.get_traceback(2, 0), [(1, 0)])


class TestScoringMatrix4D(unittest.TestCase):
    """
    Test an MSA scoring matrix for 4 sequences.
    """

    def setUp(self):
        self.sequences = ["ABC", "ABD", "ABE", "ABF"]
        self.scoring_matrix = ScoringMatrix(self.sequences)

    def test_get_score(self):
        """
        Test getting a score from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0, 0), 0)

        self.scoring_matrix.set_score(0, 0, 0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0, 0), 1)

    def test_get_max_score(self):
        """
        Test getting the maximum score from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_max_score(), 0)

        self.scoring_matrix.set_score(0, 0, 0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_max_score(), 1)

        self.scoring_matrix.set_score(1, 1, 1, 1, score=2)

        self.assertEqual(self.scoring_matrix.get_max_score(), 2)

    def test_get_max_index(self):
        """
        Test getting the maximum index from the scoring matrix.
        """
        self.assertEqual(self.scoring_matrix.get_max_index(), (0, 0, 0, 0))

        self.scoring_matrix.set_score(1, 1, 1, 1, score=2)

        self.assertEqual(self.scoring_matrix.get_max_index(), (1, 1, 1, 1))

    def test_get_traceback(self):
        """
        Test getting the traceback from the scoring matrix.
        """
        self.scoring_matrix.add_traceback(0, 0, 0, 0, traceback_direction=(1, 1, 1, 1))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0, 0), [(1, 1, 1, 1)])

        self.scoring_matrix.add_traceback(0, 0, 0, 0, traceback_direction=(2, 2, 2, 2))

        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0, 0), [(1, 1, 1, 1), (2, 2, 2, 2)])

    def test_mixed_set(self):
        """
        Test that setting both scores and tracebacks don't cause issues.
        """
        self.scoring_matrix.set_score(0, 0, 0, 0, score=1)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0, 0), [])
        self.assertEqual(self.scoring_matrix[0, 0, 0, 0], 1)
        self.assertEqual(self.scoring_matrix[0, 0, 0, 1], 0)

        self.scoring_matrix.add_traceback(0, 0, 0, 0, traceback_direction=(1, 1, 1, 1))

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0, 0), 1)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0, 0), [(1, 1, 1, 1)])
        self.assertEqual(self.scoring_matrix[0, 0, 0, 0], 1)
        self.assertEqual(self.scoring_matrix[0, 0, 0, 1], 0)

        self.scoring_matrix.set_score(0, 0, 0, 0, score=2)

        self.assertEqual(self.scoring_matrix.get_score(0, 0, 0, 0), 2)
        self.assertEqual(self.scoring_matrix.get_traceback(0, 0, 0, 0), [(1, 1, 1, 1)])
        self.assertEqual(self.scoring_matrix[0, 0, 0, 0], 2)
        self.assertEqual(self.scoring_matrix[0, 0, 0, 1], 0)
