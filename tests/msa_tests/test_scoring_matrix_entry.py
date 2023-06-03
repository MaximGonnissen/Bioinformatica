import unittest

from msa.scoring_matrix.scoring_matrix import ScoringMatrixEntry


class TestScoringMatrixEntry(unittest.TestCase):
    """
    Test ScoringMatrixEntry class.
    """

    def setUp(self) -> None:
        """
        Set up test cases.
        """
        self.entry = ScoringMatrixEntry(0, [])
        self.other_entry = ScoringMatrixEntry(1, [1, 2, 3])

    def test_equals(self) -> None:
        """
        Test equals method.
        """
        self.assertTrue(self.entry == self.entry)
        self.assertFalse(self.entry == self.other_entry)

        self.assertTrue(self.entry == (0, []))
        self.assertFalse(self.entry == (1, [1, 2, 3]))

        self.assertTrue(self.entry == 0)
        self.assertFalse(self.entry == 1)

        self.assertFalse(self.entry == None)
        self.assertFalse(self.entry == "test")

    def test_not_equals(self) -> None:
        """
        Test not_equals method.
        """
        self.assertFalse(self.entry != self.entry)
        self.assertTrue(self.entry != self.other_entry)

        self.assertFalse(self.entry != (0, []))
        self.assertTrue(self.entry != (1, [1, 2, 3]))

        self.assertFalse(self.entry != 0)
        self.assertTrue(self.entry != 1)

        self.assertTrue(self.entry != None)
        self.assertTrue(self.entry != "test")

    def test_get_item(self) -> None:
        """
        Test get_item method.
        """
        self.assertEqual(self.entry[0], 0)
        self.assertEqual(self.entry[1], [])

        self.assertEqual(self.other_entry[0], 1)
        self.assertEqual(self.other_entry[1], [1, 2, 3])

        with self.assertRaises(IndexError):
            self.entry[2]

        with self.assertRaises(IndexError):
            self.other_entry[2]

    def test_set_item(self) -> None:
        """
        Test set_item method.
        """
        self.entry[0] = 1
        self.entry[1] = [1, 2, 3]
        self.assertEqual(self.entry[0], 1)
        self.assertEqual(self.entry[1], [1, 2, 3])

        self.other_entry[0] = 0
        self.other_entry[1] = []
        self.assertEqual(self.other_entry[0], 0)
        self.assertEqual(self.other_entry[1], [])

        with self.assertRaises(IndexError):
            self.entry[2] = 1

        with self.assertRaises(IndexError):
            self.other_entry[2] = 1

    def test_str(self) -> None:
        """
        Test str method.
        """
        self.assertEqual(str(self.entry), "(0, [])")
        self.assertEqual(str(self.other_entry), "(1, [1, 2, 3])")

    def test_lt(self) -> None:
        """
        Test lt method.
        """
        self.assertTrue(self.entry < self.other_entry)
        self.assertFalse(self.other_entry < self.entry)

        self.assertTrue(self.entry < 1)
        self.assertFalse(self.other_entry < 0)

        with self.assertRaises(TypeError):
            self.entry < None

        with self.assertRaises(TypeError):
            self.entry < "test"

    def test_le(self) -> None:
        """
        Test le method.
        """
        self.assertTrue(self.entry <= self.other_entry)
        self.assertFalse(self.other_entry <= self.entry)

        self.assertTrue(self.entry <= 1)
        self.assertFalse(self.other_entry <= 0)

        with self.assertRaises(TypeError):
            self.entry <= None

        with self.assertRaises(TypeError):
            self.entry <= "test"

    def test_gt(self) -> None:
        """
        Test gt method.
        """
        self.assertTrue(self.other_entry > self.entry)
        self.assertFalse(self.entry > self.other_entry)

        self.assertTrue(self.other_entry > 0)
        self.assertFalse(self.entry > 1)

        with self.assertRaises(TypeError):
            self.entry > None

        with self.assertRaises(TypeError):
            self.entry > "test"

    def test_ge(self) -> None:
        """
        Test ge method.
        """
        self.assertTrue(self.other_entry >= self.entry)
        self.assertFalse(self.entry >= self.other_entry)

        self.assertTrue(self.other_entry >= 0)
        self.assertFalse(self.entry >= 1)

        with self.assertRaises(TypeError):
            self.entry >= None

        with self.assertRaises(TypeError):
            self.entry >= "test"

    def test_add(self) -> None:
        """
        Test add method.
        """
        self.assertEqual(self.entry + self.other_entry, 1)
        self.assertEqual(self.entry + 1, 1)
        self.assertEqual(self.other_entry + 0, 1)

        self.assertEqual(self.entry + self.entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry + 0, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry + self.other_entry, ScoringMatrixEntry(2))

        with self.assertRaises(TypeError):
            self.entry + None

        with self.assertRaises(TypeError):
            self.entry + "test"

    def test_sub(self) -> None:
        """
        Test sub method.
        """
        self.assertEqual(self.entry - self.other_entry, -1)
        self.assertEqual(self.entry - 1, -1)
        self.assertEqual(self.other_entry - 0, 1)

        self.assertEqual(self.entry - self.entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry - 0, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry - self.other_entry, ScoringMatrixEntry(0))

        with self.assertRaises(TypeError):
            self.entry - None

        with self.assertRaises(TypeError):
            self.entry - "test"

    def test_mul(self) -> None:
        """
        Test mul method.
        """
        self.assertEqual(self.entry * self.other_entry, 0)
        self.assertEqual(self.entry * 1, 0)
        self.assertEqual(self.other_entry * 0, 0)

        self.assertEqual(self.entry * self.entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry * 0, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry * self.other_entry, ScoringMatrixEntry(1))

        with self.assertRaises(TypeError):
            self.entry * None

        with self.assertRaises(TypeError):
            self.entry * "test"

    def test_truediv(self) -> None:
        """
        Test truediv method.
        """
        self.assertEqual(self.entry / self.other_entry, 0)
        self.assertEqual(self.entry / 1, 0)
        self.assertEqual(self.other_entry / 1, 1)

        self.assertEqual(self.entry / self.other_entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry / 1, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry / self.other_entry, ScoringMatrixEntry(1))

        with self.assertRaises(ZeroDivisionError):
            self.entry / 0

        with self.assertRaises(ZeroDivisionError):
            self.other_entry / 0

        with self.assertRaises(TypeError):
            self.entry / None

        with self.assertRaises(TypeError):
            self.entry / "test"

    def test_floordiv(self) -> None:
        """
        Test floordiv method.
        """
        self.assertEqual(self.entry // self.other_entry, 0)
        self.assertEqual(self.entry // 1, 0)
        self.assertEqual(self.other_entry // 1, 1)

        self.assertEqual(self.entry // self.other_entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry // 1, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry // self.other_entry, ScoringMatrixEntry(1))

        with self.assertRaises(ZeroDivisionError):
            self.entry // 0

        with self.assertRaises(ZeroDivisionError):
            self.other_entry // 0

        with self.assertRaises(TypeError):
            self.entry // None

        with self.assertRaises(TypeError):
            self.entry // "test"

    def test_mod(self) -> None:
        """
        Test mod method.
        """
        self.assertEqual(self.entry % self.other_entry, 0)
        self.assertEqual(self.entry % 1, 0)
        self.assertEqual(self.other_entry % 1, 0)

        self.assertEqual(self.entry % self.other_entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry % 1, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry % self.other_entry, ScoringMatrixEntry(0))

        with self.assertRaises(ZeroDivisionError):
            self.entry % 0

        with self.assertRaises(ZeroDivisionError):
            self.other_entry % 0

        with self.assertRaises(TypeError):
            self.entry % None

        with self.assertRaises(TypeError):
            self.entry % "test"

    def test_pow(self) -> None:
        """
        Test pow method.
        """
        self.assertEqual(self.entry ** self.other_entry, 0)
        self.assertEqual(self.entry ** 1, 0)
        self.assertEqual(self.other_entry ** 0, 1)

        self.assertEqual(self.entry ** self.other_entry, ScoringMatrixEntry(0))
        self.assertEqual(self.entry ** 1, ScoringMatrixEntry(0))
        self.assertEqual(self.other_entry ** self.other_entry, ScoringMatrixEntry(1))

        with self.assertRaises(TypeError):
            self.entry ** None

        with self.assertRaises(TypeError):
            self.entry ** "test"
