from typing import Union, Tuple, List, Any, Iterable

import numpy as np


class ScoringMatrixEntry:
    """
    Class for a scoring matrix entry.

    Primarily meant to get around an issue with numpy fill.
    """

    def __init__(self, score: Union[int, float] = 0, traceback: List[Any] = None):
        """
        Initialise the scoring matrix entry.
        :param score: Score for the entry.
        :param traceback: Traceback for the entry.
        """
        self.score = score
        self.traceback = traceback or []

    def __getitem__(self, *args) -> Any:
        """
        Get the score.
        Required to get around an issue with finding the max score.
        """
        return self.score

    def __str__(self):
        """
        Get the string representation of the scoring matrix entry.
        :return: String representation of the scoring matrix entry.
        """
        return f"({self.score}, {self.traceback})"

    def __repr__(self):
        """
        Get the representation of the scoring matrix entry.
        :return: Representation of the scoring matrix entry.
        """
        return f"ScoringMatrixEntry({self.score}, {self.traceback})"

    def __eq__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entries are equal.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entries are equal.
        """
        return self.score == other.score

    def __ne__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entries are not equal.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entries are not equal.
        """
        return self.score != other.score

    def __lt__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entry is less than the other.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entry is less than the other.
        """
        return self.score < other.score

    def __le__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entry is less than or equal to the other.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entry is less than or equal to the other.
        """
        return self.score <= other.score

    def __gt__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entry is greater than the other.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entry is greater than the other.
        """
        return self.score > other.score

    def __ge__(self, other: 'ScoringMatrixEntry') -> bool:
        """
        Check if the scoring matrix entry is greater than or equal to the other.
        :param other: Other scoring matrix entry.
        :return: Whether the scoring matrix entry is greater than or equal to the other.
        """
        return self.score >= other.score

    def __add__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Add the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Sum of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score + other.score, self.traceback + other.traceback)

    def __sub__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Subtract the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Difference of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score - other.score, self.traceback + other.traceback)

    def __mul__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Multiply the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Product of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score * other.score, self.traceback + other.traceback)

    def __truediv__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Divide the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Quotient of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score / other.score, self.traceback + other.traceback)

    def __floordiv__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Divide the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Quotient of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score // other.score, self.traceback + other.traceback)

    def __mod__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Modulo the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Modulo of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score % other.score, self.traceback + other.traceback)

    def __pow__(self, other: 'ScoringMatrixEntry') -> 'ScoringMatrixEntry':
        """
        Power the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Power of the scoring matrix entries.
        """
        return ScoringMatrixEntry(self.score ** other.score, self.traceback + other.traceback)

    def __cmp__(self, other):
        """
        Compare the scoring matrix entries.
        :param other: Other scoring matrix entry.
        :return: Comparison of the scoring matrix entries.
        """
        if self.score < other.score:
            return -1
        elif self.score > other.score:
            return 1
        else:
            return 0


class ScoringMatrix:
    """
    N-dimensional scoring matrix for MSA algorithms.
    """

    def __init__(self, sequences: List[str], *args, **kwargs):
        """
        Initialise the scoring matrix.
        :param sequences: List of sequences.
        """
        self.sequences = sequences
        self.shape = tuple(len(sequence) + 1 for sequence in sequences)
        np_from_py_wrapper = np.frompyfunc(ScoringMatrixEntry, 0, 1)
        self.matrix = np.empty(self.shape, dtype=object)
        np_from_py_wrapper(self.matrix)

    def __getitem__(self, item: Union[int, Tuple[int, ...]]) -> ScoringMatrixEntry:
        """
        Get an item from the matrix.
        :param item: Index of the item.
        :return: The item.
        """
        return self.matrix[item]

    def __setitem__(self, key: Union[int, Tuple[int, ...]], value: ScoringMatrixEntry) -> None:
        """
        Set an item in the matrix.
        :param key: Index of the item.
        :param value: Value to set the item to.
        """
        if isinstance(value, int) or isinstance(value, float):
            value = ScoringMatrixEntry(value)
        elif isinstance(value, list):
            value = ScoringMatrixEntry(traceback=value)
        elif not isinstance(value, ScoringMatrixEntry):
            raise TypeError('Value must be an integer, float, list or ScoringMatrixEntry.')
        self.matrix[key] = value

    def __iter__(self) -> Iterable:
        """
        Iterate over the matrix.
        :return: Iterator for the matrix.
        """
        return iter(self.matrix)

    def __str__(self) -> str:
        """
        Convert the matrix to a string.
        :return: String representation of the matrix.
        """
        return str(self.matrix)

    def __repr__(self) -> str:
        """
        Get the representation of the matrix.
        :return: Representation of the matrix.
        """
        return str(self.matrix)

    def get_score(self, *args) -> Union[int, float]:
        """
        Get the score for a matrix entry.
        :param args: Index of the matrix entry.
        :return: Score for the matrix entry.
        """
        return self[args].score

    def get_traceback(self, *args) -> List:
        """
        Get the traceback for a matrix entry.
        :param args: Index of the matrix entry.
        :return: Traceback for the matrix entry.
        """
        return self[args].traceback

    def set_score(self, *args, score: Union[int, float]):
        """
        Set the score for a matrix entry.
        :param args: Index of the matrix entry.
        :param score: Score for the matrix entry.
        """
        self[args].score = score

    def set_traceback(self, *args, traceback: List):
        """
        Set the traceback for a matrix entry.
        :param args: Index of the matrix entry.
        :param traceback: Traceback for the matrix entry.
        """
        self[args].traceback = traceback

    def add_traceback(self, *args, traceback_direction: Tuple[int, ...]):
        """
        Add to the traceback for a matrix entry.
        :param args: Index of the matrix entry.
        :param traceback_direction: Direction to add to the traceback, relative coordinates.
        """
        self[args].traceback.append(traceback_direction)

    def get_max_index(self) -> Tuple[int, ...]:
        """
        Get the index of the maximum score in the matrix.
        :return: Index of the maximum score in the matrix.
        """
        return np.unravel_index(np.argmax(self.matrix), self.shape)

    def max_score(self) -> Union[int, float]:
        """
        Get the maximum score in the matrix.
        :return: Maximum score in the matrix.
        """
        return np.max(self.matrix).score

    def corner_index(self) -> Tuple[int, ...]:
        """
        Get the index of the corner of the matrix.
        :return: Index of the corner of the matrix.
        """
        return tuple(dim - 1 for dim in self.shape)
